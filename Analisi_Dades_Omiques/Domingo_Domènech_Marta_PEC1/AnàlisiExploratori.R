# Anàlisi exploratori de les dades Phospho.
# El dataframe utilitzat per l'exploració prové de l'objecte creat per la classe SummarizedExperiment.
# Per tant, primer de tot importarem l'objecte "se_object":

load("se_object_Phospho.rda")


# Definició de les dades del objecte a la variable "Phospho".
Phospho <- assays(se_object)$abundance # Defining the data on a variable.
Phospho <- data.frame(Phospho)


# Explorem les dimensions i els noms de les mostres:
dim(Phospho)
names(Phospho)


# Exploració estadístiques per columna amb la funció summary convertint-lo en un dataframe per comoditat visual.
data.frame(round(apply(Phospho,2, summary))) 



# Executem una primera visualització gràfica amb un histograma amb les dades sense tractar.

# Histogrames de les dades originals
dev.new()
histo <- par(mfrow = c(3, 4))                       # Creacció de la matriu de gràfiques. Cada gràfica per cada mostra. (3 files i 4 columnes)
for (i in 1:ncol(Phospho)) {                        # Per cada columna del dataset (mostres)... 
  hist(Phospho[, i], main = names(Phospho)[i])      # ...Representa un histograma amb els valors de totes les files (variables) vs la columna (mostra) en concret. Anomena l'histograma amb el nom de la mostra.
}
par(histo)                                          # representa la matriu d'histogrames 



# Observem que les dades necessiten una transformació logaritmica. 
# La diferència entre els valors dels mínims i dels màxims del summary() també ens ho ha indicat.

Phospho_log <- Phospho
Phospho_log <- log(Phospho_log + 1) # Funció del logaritme.



# Tornem a representar els histogrames:
dev.new()
histo_log <- par(mfrow=c(3,4)) # Creacció de la matriu de gràfiques. Cada gràfica per cada mostra. (3 files i 4 columnes)

for (i in 1:ncol(Phospho_log))
  hist(Phospho_log[,i], main = names(Phospho_log)[i])

par(histo_log)



# Representació boxplot per visualitzar totes les mostres alhora.
# Boxplot amb les dades tractades pel logaritme (Phospho_log).
dev.new()
par(mar=c(7, 4, 4, 2))                                                                         # Mides de la graella del boxplot.
groupColor <- c(rep("red", 6), rep("blue", 6))                                                 # Definició dels colors segone el nombre de msotres (12) repartits entre els dos grups (MSS i PD)

boxplot(Phospho_log, col=groupColor, main="Expression values for 12 samples, two groups",      # Execució del boxplot.
        xlab="", 
        ylab="Expression", las=2, cex.axis=0.7, cex.main=0.7)
mtext(                                                                                        # Definició del disseny del títol de l'eix de les x per evitar el solapament amb els noms de les mostres.
  "Mètodes", 
  side = 1, 
  line = 5
)



# Un cop hem visualitzat les dades gràficament, procedim a centralitzar-les per poder dur a terme l'anàlisi multivariant amb PCAs.

# Centralització i normalització de les dades
Phospho_log_centered <- scale(Phospho_log, center = TRUE, scale = TRUE)

library(ggplot2)      # Llibreria per la representació de gràfiques més elaborada
library(dplyr)        # Llibreria per executar la funció prcomp


# Anàlisi multivariant: PCA --> prcomp
pca <- prcomp(t(Phospho_log_centered), scale. = TRUE)       # Degut a que hi ha més variables que mostres s'ha d'aplicar la funció t() que transposar la matriu de les dades d'entrada.

# Visualització del PCA
pca_data <- data.frame(pca$x)
pca_data$method <- col_data$method

dev.new()                                                    # Per crear una nova finestra.
ggplot(pca_data, aes(x = PC1, y = PC2, color = method)) +
  geom_point(size = 3) +
  labs(title = "PCA of phosphopeptides", x = "PC1", y = "PC2") +
  theme_minimal()




# Estudi de l'efecte batch:

# Per evitar confusions degut a l'efecte batch, és a dir evitar variacions sitemàtiques les quals no estan 
#relacionades amb les variables d'interès, sino amb amb el processament de les mostres, es duen a terme una sèrie d'anàlisis. 

# En aquest cas només es realitzaran anàlisis de detecció.

dev.new()
X <- Phospho_log_centered
# Funció plotPCA
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, cex4text=0.8)
{
  pcX <- prcomp(t(X), scale.=scale)  # Utilitzar t(X) per transposar la matriu
  loads <- round(pcX$sdev^2/sum(pcX$sdev^2)*100, 1)
  xlab <- c(paste("PC1", loads[1], "%"))
  ylab <- c(paste("PC2", loads[2], "%"))
  if (is.null(colors)) colors <- 1
  plot(pcX$x[, 1:2], xlab=xlab, ylab=ylab, col=colors,  
       xlim=c(min(pcX$x[, 1])-10, max(pcX$x[, 1])+10),
       ylim=c(min(pcX$x[, 2])-10, max(pcX$x[, 2])+10))
  text(pcX$x[, 1], pcX$x[, 2], labels, pos=3, cex=cex4text)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

# Aplicar la funció plotPCA al dataset Phospho_log_centered.
labels <- col_data$sample
colors <- as.numeric(factor(col_data$method))
dataDesc <- "Phospho_log_centered dataset"

plotPCA(Phospho_log_centered, labels=labels, colors=colors, dataDesc=dataDesc)



# Per últim executarem un estudi basat en distàncies amb la funció heatmap(). 
# Aquest estudi consisteix en calcular una matriu de distàncies i visualitzant-la mitjançant un mapa de colors.

dev.new()
manDist <- dist(t(X))
heatmap (as.matrix(manDist), col=heat.colors(16))

















