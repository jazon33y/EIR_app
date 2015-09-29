library(shiny)

cut_it_up <- function(a) {
  return(unlist(strsplit(a,split = '')))
}

plot_sequences <- function(x = NA,y = 'empty',z = 'empty') {
  len = .01
  seq_len = length(cut_it_up(x))
  iter = (1 / seq_len)
  text_display = 'on'
  border = NULL
  lwd = 2
  if (seq_len > 150) {
    border = NA;lwd = NULL;text_display = 'off'
  }
  par(mar = c(0,0,0,0))
  plot(
    c(0, 1.1), c(.3, .7), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n'
  )
  for (i in cut_it_up(x)) {
    if (i == 'A') {
      col = 'tomato2'
    }
    else if (i == 'T') {
      col = 'royalblue1'
    }
    else if (i == 'G') {
      col = 'seagreen2'
    }
    else if (i == 'C') {
      col = 'khaki2'
    }
    else {
      col = 'grey'
    }
    len = len + iter
    polygon(
      x = c(
        len - (.06 / (seq_len / 10)),len - (.05 / (seq_len / 10)),
        len - (.06 / (seq_len / 10)),len + (.04 / (seq_len / 10)),
        len + (.05 / (seq_len / 10)),len + (.04 / (seq_len / 10)),
        len - (.06 / (seq_len / 10))
      ),
      y = c(0.41,0.5,0.59,0.59,0.5,0.41,0.41),col = col,lwd = lwd,border =
        border
    )
    if (text_display == 'on') {
      text(
        x = len, y = .5, paste(i),col = 'black',cex = 2.5 / (seq_len / 10)
      )
    }
  }
  if ((x != 'empty') && (y != 'empty')) {
    len = 0.01
    polygon(
      x = c(((iter * y) + len) - (.06 / (seq_len / 10)),
            ((iter * y) + len) - (.06 / (seq_len / 10)),
            ((iter * z) + len) + (.04 / (seq_len / 10)),
            ((iter * z) + len) + (.04 / (seq_len / 10))
      ),
      y = c(.38,.62,.62,.38),col = rgb(.2, .2, .2, .2),border = border,lty =
        2
    )
  }
}


eir <- function(chrom,start_loc,indel) {
  library(BSgenome.Hsapiens.UCSC.hg19)
  getSeq_char <- function(genome,chrom = 'chrM',loc1,loc2) {
    seq = as.character(getSeq(genome,chrom,loc1,loc2))
    seq = strsplit(seq,split = '')
    return(unlist(seq))
  }
  
  genome <- BSgenome.Hsapiens.UCSC.hg19
  ref = getSeq_char(genome,paste0('chr',chrom),start_loc - 1000,start_loc + 1000)
  
  s = ref
  m = length(indel)
  d = u = 1001#start_loc
  k = 1
  while (indel[k] == s[d + m]) {
    #   cat("k = ",k,"d + m = ", d+m,"\nindel[k] = ",indel[k],"\t s[d+m] = ",s[d+m],"\n")
    d = d + 1
    k = k + 1
    if (k == (m + 1)) {
      k = 1
    }
    #   cat("k = ",k,"d + m = ", d+m,"\nindel[k] = ",indel[k],"\t s[d+m] = ",s[d+m],"\n\n")
  }
  k = m
  while (indel[k] == s[u - 1]) {
    u = u - 1
    k = k - 1
    if (k < 1) {
      k = m
    }
  }
  u = 1001 - u
  d = d - 1001
  return(c((start_loc - u),start_loc + d))
  #   return(c(u,d))
}

indel_data = read.delim(
  "chr22_parsed",header = FALSE
)
colnames(indel_data) <-
  c('Chr','Left','Right','RS_id','Freq','Allele')

library(BSgenome.Hsapiens.UCSC.hg19)
getSeq_char <- function(genome,chrom = 'chrM',loc1,loc2) {
  seq = as.character(getSeq(genome,chrom,loc1,loc2))
  seq = strsplit(seq,split = '')
  return(unlist(seq))
}

genome <- BSgenome.Hsapiens.UCSC.hg19

shinyServer(function(input, output) {
  output$plot <- renderPlot({
    m = input$chr
    loc = as.integer(input$loc)
    indel = input$allele
    indel = cut_it_up(as.character(indel))
    answer = eir(m,loc,indel)
    #     print(m)
    #     print(loc)
    #     print(answer)
    #     print(indel)
    
    l = answer[1]
    r = answer[2]
    offset = 10
    seql = l - offset
    seqr = r + offset
    
    # size = r - l
    # print( size + 1 )
    
    get_seq = getSeq_char(genome,paste0('chr',m),seql,seqr)
    plot_sequences(get_seq,offset + 1,length(cut_it_up(get_seq)) - offset)
    
  })
  
  # Filter data based on selections
  output$table <- renderDataTable({
    data <- indel_data # eir table
    CHROM <- paste0(c('chr',as.character(input$chr)),collapse = '')
    data <- data[data$Chr == CHROM,]
    DL = as.integer(data$Left)
    DR = as.integer(data$Right)
    Q = as.integer(as.character(input$loc))
    data <- data[(Q >= DL) & (Q <= DR),]
    data <- data[data$Allele == as.character(input$allele),]
    data
  })
  
})
