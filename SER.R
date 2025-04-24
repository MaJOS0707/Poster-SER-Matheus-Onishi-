#-------------------------------------------------------------------------------
# Bliblioteca
library(BayesLCA)
library(ggplot2)
library(dplyr)  
library(tidyr) 
library(coda)
library(tibble)
library(spdep)
library(tmap)
library(sf)
library(readr)
library(leaflet)
library(leafsync)
library(htmlwidgets)
library(magrittr)

library(spdep)

library(RColorBrewer)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Itens salvos


baseCisp_binaria_20a23 = read_csv2("baseCisp_binaria_20a23.csv")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# LCA
LCA = blca.gibbs(baseCisp_binaria_20a23[,-1],2,burn.in = 2000, iter = 12000)


# Probabilidade das classes
LCA$classprob

# Probabilidade dos itens
LCA$itemprob

# Pegando as Cadeias 
mcmc_class <- as.mcmc(LCA$samples$classprob)
mcmc_itens_G1 = as.mcmc(LCA$samples$itemprob[,1,])
mcmc_itens_G2 = as.mcmc(LCA$samples$itemprob[,2,])


# Calcular os intervalos de credibilidade de 95%
intervals_class <- HPDinterval(mcmc_class, prob = 0.95)
intervals_itens_G1 =  HPDinterval(mcmc_itens_G1, prob = 0.95)
intervals_itens_G2=  HPDinterval(mcmc_itens_G2, prob = 0.95)

# Exibir os intervalos
intervals_class
intervals_itens_G1
intervals_itens_G2

class_posteriori =  function(base,LCA,G){
  # Número de classes e itens
  num_classes <- G
  
  # Extraindo probabilidades dos itens e proporções de classe
  class_prior <- LCA$classprob
  item_probs <- LCA$itemprob
  # Inicializando matriz para armazenar as probabilidades posteriores
  posterior_probs_matrix <- matrix(0, nrow = nrow(base), ncol = num_classes)
  # Loop sobre as observações
  for (i in 1:nrow(base)) {
    response_i <- as.numeric(base[i, -c(1,2)])  # Obter respostas da observação i
    
    # Calcular as probabilidades condicionais para cada classe
    class_cond_probs <- numeric(num_classes)
    for (k in 1:num_classes) {
      p_y_given_class_k <- prod(response_i * item_probs[k, ] + (1 - response_i) * (1 - item_probs[k, ]))
      class_cond_probs[k] <- p_y_given_class_k * class_prior[k]
    }
    
    # Normalizar para obter P(Classe = k | Y = Y_i) e armazenar na matriz
    posterior_probs_matrix[i, ] <- class_cond_probs / sum(class_cond_probs)
  }
  # A matriz posterior_probs_matrix contém uma linha por observação e uma coluna por classe.
  return(data.frame(Prob_G1 = posterior_probs_matrix[,1],
                    Prob_G2 = posterior_probs_matrix[,2]))
}

probs_Class_posterior = class_posteriori(baseCisp_binaria_20a23,LCA,G=2)


# Juntando as bases
baseCisp_binaria_20a23 = cbind(baseCisp_binaria_20a23,probs_Class_posterior)

#----------------------------------------------------------------------------------------------
# Mapa

# Carregar o shapefile do RJ
rj = st_read(dsn = "lm_cisp_bd.shp")
bbox <- st_bbox(mapa)  # mapa = seu objeto sf com geometria
mapa <- rj |> 
  left_join(baseCisp_binaria_20a23, by = "cisp")

# Criando uma paleta de cores
pal <- colorNumeric("YlOrRd", domain = 0:1)

mapa_crime <- leaflet(mapa, options = leafletOptions(zoomControl = FALSE)) |> 
  addTiles() |> 
  addPolygons(
    fillColor = ~pal(Prob_G1),
    weight = 1,
    color = "white",
    fillOpacity = 0.8,
    label = ~paste0("CISP: ", cisp, "<br>Prob. Classe 1: ", round(Prob_G1, 5)),
    highlightOptions = highlightOptions(weight = 2, color = "black", bringToFront = TRUE)
  ) |> 
  addLegend(
    pal = pal, 
    values = ~Prob_G1, 
    opacity = 0.7, 
    title = "Probabilidade Classe 1",
    position = "bottomright"     # <- legenda inferior direita
  ) |> 
  addControl(
    html = "<h3 style='color:black;text-align:center;'></h3>", 
    position = "topright"
  )
mapa_crime





#############################################################

# Certificar geometria válida
mapa <- st_make_valid(mapa)

# Converter para objeto Spatial
mapa_sp <- as(mapa, "Spatial")

# Definir vizinhança e matriz de pesos
viz <- poly2nb(mapa_sp)
lw <- nb2listw(viz, style = "W", zero.policy = TRUE)

# Teste de Moran Local (LISA)
lisa <- localmoran(mapa$Prob_G1, lw, zero.policy = TRUE)

# Adicionar estatísticas ao shapefile
mapa$I_local <- round(lisa[, "Ii"], 3)
mapa$p_value <- lisa[, "Pr(z != E(Ii))"]

# Definir média para classificação LISA
media_prob <- mean(mapa$Prob_G1)

# Classificação LISA
mapa$LISA_Cat <- "Não Significativo"
mapa$LISA_Cat[mapa$Prob_G1 > media_prob & mapa$I_local > 0 & mapa$p_value < 0.05] <- "High-High"
mapa$LISA_Cat[mapa$Prob_G1 < media_prob & mapa$I_local > 0 & mapa$p_value < 0.05] <- "Low-Low"
mapa$LISA_Cat[mapa$Prob_G1 > media_prob & mapa$I_local < 0 & mapa$p_value < 0.05] <- "High-Low"
mapa$LISA_Cat[mapa$Prob_G1 < media_prob & mapa$I_local < 0 & mapa$p_value < 0.05] <- "Low-High"

# Converter para fator com ordem personalizada
mapa$LISA_Cat <- factor(mapa$LISA_Cat, 
                        levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Não Significativo"))

# Paleta de cores manual
cores_lisa <- c("High-High" = "darkred",
                "Low-Low" = "darkblue",
                "High-Low" = "orange",
                "Low-High" = "lightblue",
                "Não Significativo" = "grey")

pal_lisa <- colorFactor(palette = cores_lisa, domain = levels(mapa$LISA_Cat))

# Criar mapa interativo com leaflet
mapa_lisa <- leaflet(mapa) %>%
  addTiles() %>%
  addPolygons(
    fillColor = ~pal_lisa(LISA_Cat),
    weight = 1,
    color = "white",
    fillOpacity = 0.8,
    label = ~paste0("CISP: ", cisp, 
                    "<br>Categoria LISA: ", LISA_Cat,
                    "<br>I Local: ", I_local,
                    "<br>p-valor: ", round(p_value, 3)),
    highlightOptions = highlightOptions(weight = 2, color = "black", bringToFront = TRUE)
  )|> 
  addLegend(
    pal =  pal_lisa, 
    values = ~levels(mapa$LISA_Cat), 
    opacity = 0.7, 
    title = "LISA - Moran Local",
    position = "bottomright"     # <- legenda inferior direita
  ) |> 
  addControl(
    html = "<h3 style='color:black;text-align:center;'></h3>", 
    position = "topright"
  )


# Visualizar
mapa_lisa
