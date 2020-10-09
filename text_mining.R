library(udpipe)
library(tidyverse)
library(tm)
library(stargazer)
library(qdap)
library(lattice)
library(igraph)
library(ggraph)
library(gplots)
library(ggplot2)
library(tidygraph)
library(ggpubr)
library(factoextra)
library(wordcloud)
library(quanteda)
library(corrplot)

setwd("~/R/...")

gp <- read_csv("keywords.csv",
               col_types = cols(date = col_date(format = "%Y"),
                                keywords = col_skip(), 
                                similarity = col_skip()))

corpus <- SimpleCorpus(VectorSource(gp$text))


#special chars
removeSpecialChars <- function(x) gsub("[^a-zA-Z0-9 ]","", x)
corpus <- tm_map(dfCorpus, removeSpecialChars)


#stem
corpus_stem <- tm_map(corpus, stemDocument)


#dtm
dtm <- DocumentTermMatrix(corpus_stem)
dtm <- removeSparseTerms(dtm, .97)
inspect(dtm)


#inspect most frequent words
findFreqTerms(dtm, lowfreq=500)
#[1] "access"  "develop" "project" "region"  "sustain" "tourism" "area"   
#[8] "local"   "product" "promot"  "citi"    "tourist" "visitor" "cultur" 
#[15] "heritag"

freq_words <- colnames(t(findMostFreqTerms(dtm, n = 20, 
                                           INDEX = rep(1, DTM$nrow))[[1]]))


#create frequency table
m <- tm_map(corpus_stem, removeWords, c("tourism", "tourist", "visitor"))
m <- DocumentTermMatrix(m)
m <- removeSparseTerms(m, 0.999)
m <- t(as.matrix(m))

freq_table <- data.frame(term = rownames(m), 
                         freq = rowSums(m), 
                         row.names = NULL)

freq_table <- freq_table[order(-freq_table$freq),][1:20,]
freq_table

#        term freq
#29   develop  989
#79    region  887
#75   project  806
#165    local  780
#496   cultur  761
#93   sustain  721
#229     citi  615
#115     area  576
#185   promot  545
#182  product  541
#1     access  509
#535  heritag  509
#3      activ  497
#1057   natur  461
#65       new  455
#166    manag  444
#86    servic  440
#132   destin  414
#175    offer  408
#102      use  394


#frequency plot
freq_plot <- ggplot(freq_table, aes(x = reorder(term, -freq), freq)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(x = "Terms", y = "Frequency", title = "Frequent terms") +
  geom_text(aes(label = freq), vjust = -0.5, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

freq_plot


#co-occurence heatmap
dtm_freq_words <- dtm[, Terms(dtm) %in% freq_words]
dtm_freq_words

m_freq_words <- as.matrix(dtm_freq_words)
heatmap_data <- t(m_freq_words) %*% m_freq_words
heatmap_data

diag(heatmap_data) <- 0
heatmap_data[lower.tri(heatmap_data)] <- NA

heatmap <- heatmap.2(m,
           dendrogram = "none", Colv = FALSE, Rowv = FALSE,
           scale = "none", col = brewer.pal(5, "Blues"),
           key = TRUE, density.info = "none", key.title = NA, key.xlab = "Frequency",
           trace = "none",
           main = "Term co-occurrence",
           xlab = "Term",
           ylab = "Term")


#term correlation
cor_data <- cor(m_freq_words)

cor_plot <- corrplot(cor_data, method = "square", 
                     type = "upper", tl.col = "black", 
                     order = "hclust", col = brewer.pal(n = 5, name = "RdYlBu"))


#word network
set.seed(12345)
remove = c("tourist", "tourism", "visitor")

text <- gp$text

tokens <- text %>%
  tokens(remove_punct = TRUE) %>%
  tokens_tolower() %>%
  tokens_remove(pattern = stopwords("english"), padding = FALSE)

fcmat <- fcm(toks2, context = "window", tri = FALSE)

feat <- names(topfeatures(fcmat, 30))

fcm <- fcm_select(fcmat, pattern = feat)

word_network <- textplot_network(fcm, 
                                 vertex_labelsize = 3*rowSums(fcm)/min(rowSums(fcm)),
                                 min_freq = 0.6)


## Clustering
# Building the feature matrices
dtm_tfidf <- tm::weightTfIdf(dtm)
tfidf_matrix <- as.matrix(dtm_tfidf)

#cosine distance matrix
dist_matrix <- proxy::dist(tfidf_matrix, method = "cosine")

#k-means
k <- 3

k_cluster <- kmeans(tfidf_matrix, k, nstart = 25)

k_plot <- fviz_cluster(k_cluster, tfidf_matrix, 
                       ellipse = FALSE, 
                       geom = "point")

#hierarchical clusters
h_cluster <- hcut(dist_matrix, k)

h_plot <- fviz_cluster(h_cluster, dist_matrix, 
                       ellipse = FALSE, 
                       geom = "point")

h_dendogram <- fviz_dend(h_cluster, rect = TRUE)

#test optimal value for k
data <- tfidf_matrix

elbow_h <- fviz_nbclust(data, hcut, method = "wss", k.max = 10) +
  labs(subtitle = "Elbow method")

sil <- fviz_nbclust(data, hcut, method = "silhouette", k.max = 10) +
  labs(subtitle = "Silhouette method")

#add cluster data to dataset
gp$cluster <- h_cluster$cluster

cluster_1 <- gp %>%
  select(id, source, title, type_ecst, institution, text, summary_clean, 
         location, keywords_reviewed, cluster) %>%
  filter(cluster == 1)

cluster_2 <- gp %>%
  select(id, source, title, type_ecst, institution, text, summary_clean, 
         location, keywords_reviewed, cluster) %>%
  filter(cluster == 2)

cluster_3 <- gp %>%
  select(id, source, title, type_ecst, institution, text, summary_clean, 
         location, keywords_reviewed, cluster) %>%
  filter(cluster == 3)


#cluster 1

corpus_1 <- SimpleCorpus(VectorSource(cluster_1$summary_clean))
corpus_1 <- tm_map(corpus_1, 
                   removeWords, c("tourism", "tourist", "visitor"))

#document term matrix
dtm_1 <- DocumentTermMatrix(corpus_1)
dtm_1 <- removeSparseTerms(dtm_1, .97)
inspect(dtm_1)

#wordcloud
matrix_1 <- as.matrix(dtm_1)
words_1 <- sort(colSums(matrix_1),decreasing=TRUE)
df_1 <- data.frame(word = names(words_1),freq=words_1)
set.seed(1234) # for reproducibility
wordcloud(words = df_1$word, freq = df_1$freq, min.freq = 1, 
          scale=c(4,0.5), max.words=100, random.order=FALSE, 
          rot.per=0.35, colors=brewer.pal(3, "Paired"))

#inspect most popular words
findFreqTerms(dtm_1, lowfreq=100)

mat1 <- t(matrix_1)

freq_tab1 <- data.frame(term = rownames(mat1), 
                        freq = rowSums(mat1), 
                        row.names = NULL)

freq_tab1 <- freq_tab1[order(-freq_tab1$freq),][1:10,]
freq_tab1$cluster = 1

freq_plot1 <- ggplot(freq_tab1, aes(x = reorder(term, -freq), freq)) +
  geom_bar(stat = "identity", fill = "#F8766D") +
  labs(x = "Terms", y = "Frequency", title = "Frequent terms in cluster #1") +
  geom_text(aes(label = freq), vjust = -0.5, size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#cluster 2

corpus_2 <- SimpleCorpus(VectorSource(cluster_2$summary_clean))
corpus_2 <- tm_map(corpus_2, 
                   removeWords, c("tourism", "tourist", "visitor", "tourists", "visitors"))

#document term matrix
dtm_2 <- DocumentTermMatrix(corpus_2)
dtm_2 <- removeSparseTerms(dtm_2, .97)
inspect(dtm_2)

#wordcloud
matrix_2 <- as.matrix(dtm_2)
words_2 <- sort(colSums(matrix_2),decreasing=TRUE)
df_2 <- data.frame(word = names(words_2),freq=words_2)
set.seed(1234) # for reproducibility
wordcloud(words = df_2$word, freq = df_2$freq, min.freq = 1, 
          scale=c(4,0.5), max.words=100, random.order=FALSE, 
          rot.per=0.35, colors=brewer.pal(3, "Paired"))

#inspect most popular words
findFreqTerms(dtm_2, lowfreq=50)

mat2 <- t(matrix_2)

freq_tab2 <- data.frame(term = rownames(mat2), 
                        freq = rowSums(mat2), 
                        row.names = NULL)

freq_tab2 <- freq_tab2[order(-freq_tab2$freq),][1:10,]
freq_tab2$cluster = 2

freq_plot2 <- ggplot(freq_tab2, aes(x = reorder(term, -freq), freq)) +
  geom_bar(stat = "identity", fill = "#7CAE00") +
  labs(x = "Terms", y = "Frequency", title = "Frequent terms in cluster #2") +
  geom_text(aes(label = freq), vjust = -0.5, size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#cluster 3
corpus_3 <- SimpleCorpus(VectorSource(cluster_3$summary_clean))
corpus_3 <- tm_map(corpus_3, removeWords, c("tourism", "tourist", "visitor"))

#document term matrix
dtm_3 <- DocumentTermMatrix(corpus_3)
dtm_3 <- removeSparseTerms(dtm_3, .97)
inspect(dtm_3)

#wordcloud
matrix_3 <- as.matrix(dtm_3)
words_3 <- sort(colSums(matrix_3),decreasing=TRUE)
df_3 <- data.frame(word = names(words_3),freq=words_3)
set.seed(1234) # for reproducibility
wordcloud(words = df_3$word, freq = df_3$freq, min.freq = 1, 
          scale=c(4,0.5), max.words=100, random.order=FALSE, 
          rot.per=0.35, colors=brewer.pal(3, "Paired"))

#inspect most popular words
findFreqTerms(dtm_3, lowfreq=100)

mat3 <- t(matrix_3)

freq_tab3 <- data.frame(term = rownames(mat3), 
                        freq = rowSums(mat3), 
                        row.names = NULL)

freq_tab3 <- freq_tab3[order(-freq_tab3$freq),][1:10,]
freq_tab3$cluster = 3

freq_plot3 <- ggplot(freq_tab3, aes(x = reorder(term, -freq), freq)) +
  geom_bar(stat = "identity", fill = "#00BFC4") +
  labs(x = "Terms", y = "Frequency", title = "Frequent terms in cluster #3") +
  geom_text(aes(label = freq), vjust = -0.5, size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#RAKE
#full text rake
ud_model <- udpipe_download_model(language = "english")
ud_model <- udpipe_load_model(ud_model$file_model)
x <- udpipe_annotate(ud_model, x = gp$text)
x <- as.data.frame(x)


rake <- keywords_rake(x = x, 
                      term = "token", 
                      group = c("doc_id", "paragraph_id", "sentence_id"),
                      relevant = x$upos %in% c("NOUN", "ADJ"),
                      ngram_max = 2)

rake_kw <- subset(rake, freq > 3)
rake_kw[1:20,]

#rake cluster 1
x1 <- udpipe_annotate(ud_model, x = cluster_1$text)
x1 <- as.data.frame(x1)
rake1 <- keywords_rake(x = x1, 
                       term = "token", group = c("doc_id", "paragraph_id", "sentence_id"),
                       relevant = x1$upos %in% c("NOUN", "ADJ"),
                       ngram_max = 2)

rake_kw1 <- subset(rake1, freq > 2)
rake_kw1[1:20,]

#rake cluster 2
x2 <- udpipe_annotate(ud_model, x = cluster_2$text)
x2 <- as.data.frame(x2)
rake2 <- keywords_rake(x = x2, 
                       term = "token", group = c("doc_id", "paragraph_id", "sentence_id"),
                       relevant = x2$upos %in% c("NOUN", "ADJ"),
                       ngram_max = 2)

rake_kw2 <- subset(rake2, freq > 2)
rake_kw2[1:10,]
stargazer(rake_kw2[1:20,], summary = FALSE)

#rake cluster 3
x3 <- udpipe_annotate(ud_model, x = cluster_3$text)
x3 <- as.data.frame(x3)
rake3 <- keywords_rake(x = x3, 
                       term = "token", group = c("doc_id", "paragraph_id", "sentence_id"),
                       relevant = x3$upos %in% c("NOUN", "ADJ"),
                       ngram_max = 2)

rake_kw3 <- subset(rake3, freq > 2)
rake_kw3[1:10,]
