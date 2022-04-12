
library(shiny)
library('Biostrings')
library('dplyr')
library('ggplot2')
library('seqTools')
# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel('Sequence Statistics'),

     
    sidebarLayout(
        sidebarPanel(
            fileInput("file", "Please upload a .fasta file",
                      accept = ".fasta"),
            selectInput("seq", "Sequence ID",choices='Please upload a file', FALSE),
            sliderInput("kmer",
                        "K-mer Size",
                        min = 3,
                        max = 10,
                        value = 10),
            selectInput("translation", "Translation",choices=c('Forwards, frame 1','Forwards, frame 2','Forwards, frame 3','Reverse, frame 1','Reverse, frame 2','Reverse, frame 3'), FALSE),
        
            ),

        
        mainPanel(
           
           textOutput('seq_no'),
           textOutput('Txt'),
           textOutput('GC_content'),
           textOutput('Subseq'),
           tableOutput('topKmer'),
           verbatimTextOutput('Translation'),
           plotOutput("distPlot")
        )
    )
)


server <- function(input, output, session) {
    
    inFile <- reactive({
        if (is.null(input$file)) {
            return(NULL)
        } else {
            input$file
        }
    })
    
    fasta <- reactive({
        readDNAStringSet(inFile()$datapath)
        
    })
    
    seq_name <- reactive({
        if (is.null(inFile())) {
            return(NULL)
        } else {
            names(readDNAStringSet(inFile()$datapath))
        }
    })
    
    sequence <- reactive({
        if (is.null(inFile())) {
            return(NULL)
        } else {
            paste(fasta())
        }
    })
    
    df <- reactive({
        if (is.null(inFile())) {
            return(NULL)
        } else {
            data.frame(seq_name(), sequence())
        }
    })
    
    observe({
        updateSelectInput(
            session,
            "seq",
            choices=seq_name())
        
    })

    output$seq_no<- renderText({
        if (is.null(inFile())) {
            return(NULL)
        } else {
            paste("Number of sequences in file:",nrow(df()))
        }
    })
    
    output$GC_content<- renderText({
        if (is.null(inFile())) {
            return(NULL)
        } else {
        bases <- c(strsplit(df()$sequence[df()$seq_name==input$seq], "")[[1]])
        len<-nchar(df()$sequence[df()$seq_name==input$seq])
        freq=data.frame(table(bases))
        paste("GC content:",round((( freq$Freq[freq$bases=='G'] + freq$Freq[freq$bases=='C'])/len),digits=3))
        }
        })
    
    output$Txt<- renderText({
        if (is.null(inFile())) {
            return(NULL)
        } else {
        paste("Sequence length:",nchar(df()$sequence[df()$seq_name==input$seq]))
        }
        })
    
    
    output$Subseq<- renderText({
        if (is.null(inFile())) {
            return(NULL)
        } else {
            bases2 <- c(df()$sequence[df()$seq_name==input$seq][[1]])
            n<-input$kmer
            kmer_freq<-countDnaKmers(bases2, n, 1, nchar(bases2)-(n-1))
            kmer_freq_df<-as.data.frame(kmer_freq)
            colnames(kmer_freq_df)<-'freq'
            
            G_C<-0
            kmer_freq_df_nozero <- filter(kmer_freq_df, freq > 0)
            top_5<-head(kmer_freq_df_nozero %>% arrange(desc(freq)),n=5)
            for (row in rownames(kmer_freq_df_nozero)){
                counter<-0
                for (i in strsplit(row,'')[[1]]){
                    if (i=='G' || i=='C'){
                        counter<-counter+1
                    }
                    
                }
                if (counter>G_C){
                    G_C<-counter
                    winner<-row
                }
            }
            paste("The k-mer with greatest GC content is:",winner)
            
        }
        
    })
    
    output$topKmer<- renderTable({
        if (is.null(inFile())) {
            return(NULL)
        } else {
        bases2 <- c(df()$sequence[df()$seq_name==input$seq][[1]])
        n<-input$kmer
        kmer_freq<-countDnaKmers(bases2, n, 1, nchar(bases2)-(n-1))
        kmer_freq_df<-as.data.frame(kmer_freq)
        colnames(kmer_freq_df)<-'freq'
        
        G_C<-0
        kmer_freq_df_nozero <- filter(kmer_freq_df, freq > 0)
        top_5<-head(kmer_freq_df_nozero %>% arrange(desc(freq)),n=5)
        top_5$kmer<-rownames(top_5)
        print('Top 5 k-mers:')
        top_5
        }
    }, caption="Most frequent k-mers in the sequence")
    
    output$Translation<- renderText({
        if (is.null(inFile())) {
            return(NULL)
        } else {
        fasta_all_frames <- lapply(1:3, function(pos) subseq(c(fasta(), reverseComplement(fasta())), start=pos))
        translations <-lapply(fasta_all_frames, translate,if.fuzzy.codon = 'X')
        
        if (input$translation=='Forwards, frame 1'){
            frame<-1
            direction<-1
        }
        
        else if (input$translation=='Forwards, frame 2'){
            frame<-2
            direction<-1
        }
        else if (input$translation=='Forwards, frame 3'){
            frame<-3
            direction<-1
        }
        else if (input$translation=='Reverse, frame 1'){
            frame<-1
            direction<-2
        }
        else if (input$translation=='Reverse, frame 2'){
            frame<-2
            direction<-2
        }
        else if (input$translation=='Reverse, frame 3'){
            frame<-3
            direction<-2
        }
        aa2df <- function(aa) data.frame(width=width(aa), seq=as.character(aa), names=names(aa))
        current_frame <- aa2df(translations[[frame]])
        current_translation <- current_frame$seq[current_frame$names==input$seq][direction]
        paste("Translation:",current_translation,collapse = "\n")
        }
    })
    
    output$distPlot <- renderPlot({
        if (is.null(inFile())) {
            return(NULL)
        } else {
        bases <- c(strsplit(df()$sequence[df()$seq_name==input$seq], "")[[1]])
        len<-nchar(df()$sequence[df()$seq_name==input$seq])
        freq=data.frame(table(bases))
        freq$Proportion=freq$Freq/len
        freq %>% arrange(bases)
        colours<-c('blue','black','red','green','yellow','orange')
        ggplot(freq, aes(x=bases, y=Proportion,fill = bases)) +  
            geom_bar(stat='identity' ) +
            scale_fill_manual(values = colours ) +
            theme(legend.position="none") +
            theme(text = element_text(size = 20))
        
    }
        }, height = 600, width = 500)
}

# Run the application 
shinyApp(ui = ui, server = server)
