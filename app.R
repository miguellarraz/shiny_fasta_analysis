#Libraries
library(shiny)
library('Biostrings')
library('dplyr')
library('ggplot2')
library('seqTools')
# ui
ui <- fluidPage(
    #title of tab
    titlePanel('Sequence Statistics'),

     
    sidebarLayout(
        sidebarPanel(
            #4 options, file input box, a seq id menu, a kmer size slider bar and a translation type menu
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
           #All the outputs: Number of sequences, seq length, GCcontent, kmer, top kmers table and the base content plot
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

#Server with session to make it reactive
server <- function(input, output, session) {
    #Returns NULL, so that no errors are returned before a file is uploaded
    inFile <- reactive({
        if (is.null(input$file)) {
            return(NULL)
        } else {
            input$file
        }
    })
    #Read the input, as fasta file
    fasta <- reactive({
        readDNAStringSet(inFile()$datapath)
        
    })
    #get seq names, so that menu can be updated
    seq_name <- reactive({
        if (is.null(inFile())) {
            return(NULL)
        } else {
            names(readDNAStringSet(inFile()$datapath))
        }
    })
    #Reactive, making it easier to update seq choice
    sequence <- reactive({
        if (is.null(inFile())) {
            return(NULL)
        } else {
            paste(fasta())
        }
    })
    #Df with the names and sequences in the file
    df <- reactive({
        if (is.null(inFile())) {
            return(NULL)
        } else {
            data.frame(seq_name(), sequence())
        }
    })
    #Update choices
    observe({
        updateSelectInput(
            session,
            "seq",
            choices=seq_name())
        
    })
#First output with number of rows in df which is the same as number of sequences
    output$seq_no<- renderText({
        if (is.null(inFile())) {
            return(NULL)
        } else {
            paste("Number of sequences in file:",nrow(df()))
        }
    })
    #GC content
    output$GC_content<- renderText({
        if (is.null(inFile())) {
            return(NULL)
        } else {
        #Make a list with all the bases in the sequence    
        bases <- c(strsplit(df()$sequence[df()$seq_name==input$seq], "")[[1]])
        #seq length
        len<-nchar(df()$sequence[df()$seq_name==input$seq])
        #table with all the base frqiencies
        freq=data.frame(table(bases))
        #Gc content
        paste("GC content:",round((( freq$Freq[freq$bases=='G'] + freq$Freq[freq$bases=='C'])/len),digits=3))
        }
        })
    #seq length with number of characters in sequnece
    output$Txt<- renderText({
        if (is.null(inFile())) {
            return(NULL)
        } else {
        paste("Sequence length:",nchar(df()$sequence[df()$seq_name==input$seq]))
        }
        })
    
    #get all kmers
    output$Subseq<- renderText({
        if (is.null(inFile())) {
            return(NULL)
        } else {
            #string with the bases
            bases2 <- c(df()$sequence[df()$seq_name==input$seq][[1]])
            #chosen kmer size
            n<-input$kmer
            #use count kmers to make all possible kmers and tabulate frequencies
            #Useful because I can use it for top 5 kmers as well as the greatest gc content
            kmer_freq<-countDnaKmers(bases2, n, 1, nchar(bases2)-(n-1))
            kmer_freq_df<-as.data.frame(kmer_freq)
            colnames(kmer_freq_df)<-'freq'
            #loop through all the kmers that are present at least once, split them and if its gc content is greater than the top one before, update the GC counter and record that kmer as the current winner
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
    #Top 5 kmers, simmilar method as before, but output as table
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
        #arange in descending order
        top_5<-head(kmer_freq_df_nozero %>% arrange(desc(freq)),n=5)
        top_5$kmer<-rownames(top_5)
        print('Top 5 k-mers:')
        top_5
        }
    }, caption="Most frequent k-mers in the sequence")
    
    #list with all possible transaltions of the shosen sequence, and depending on user input, one is selected
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
        #Function converting the resulting values into a df, easier to extract values from
        aa2df <- function(aa) data.frame(width=width(aa), seq=as.character(aa), names=names(aa))
        current_frame <- aa2df(translations[[frame]])
        current_translation <- current_frame$seq[current_frame$names==input$seq][direction]
        paste("Translation:",current_translation,collapse = "\n")
        }
    })
    #Base content plot
    output$distPlot <- renderPlot({
        if (is.null(inFile())) {
            return(NULL)
        } else {
        bases <- c(strsplit(df()$sequence[df()$seq_name==input$seq], "")[[1]])
        len<-nchar(df()$sequence[df()$seq_name==input$seq])
        freq=data.frame(table(bases))
        freq$Proportion=freq$Freq/len
        #arrange in alphabetical order so that bases tend to be plotted with the same colour (unless there are fuzzy bases in the data)
        freq %>% arrange(bases)
        #More than 4 colours to account for fuzzy bases
        colours<-c('blue','black','red','green','yellow','orange')
        ggplot(freq, aes(x=bases, y=Proportion,fill = bases)) +  
            geom_bar(stat='identity',color="white" , width=0.75,size=2,alpha=0.85) +
            scale_fill_manual(values = colours ) +
            geom_text(aes(label=round(Proportion,digits=3)), vjust=1.6, color="white", size=7)+
            theme(legend.position="none") +
            theme(text = element_text(size = 20)) 
        
    }
        }, height = 600, width = 500)
}

# Run the app 
shinyApp(ui = ui, server = server)
