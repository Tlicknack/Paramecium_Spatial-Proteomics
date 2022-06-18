makeWordCloud = function(df, cw){
  xxx = Corpus(VectorSource(df$Description))
  xxx <- tm_map(xxx, removePunctuation) # Remove punctuations
  xxx <- tm_map(xxx, removeWords, cw) 
  dtm <- TermDocumentMatrix(xxx)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  set.seed(42)
  wordcloud(words = d$word, freq = d$freq, min.freq = 1,
            max.words=200, random.order=FALSE, rot.per=0.35, 
            colors=brewer.pal(8, "Dark2"))
  
  # Give a dataframe with gene description and vector of "Common words" to remove
}
