# =====================================================
# CLUSTERING JERARQUICO (Ward) + ACP  -  dataset pottery
# =====================================================

# ---- 0. PAQUETES --------------------------------------
paquetes <- c("HSAUR3", "FactoMineR", "factoextra", "ggplot2")
nuevos <- paquetes[!paquetes %in% rownames(installed.packages())]
if (length(nuevos) > 0) install.packages(nuevos)
invisible(lapply(paquetes, library, character.only = TRUE))

# ---- 1. ESTILO PARA DIAPOSITIVAS ----------------------
# Colores oficiales UNAL: verde, rojo, azul oscuro, amarillo, azul claro, verde claro
col_unal <- c("#94B43B", "#A61C31", "#35408F", "#F3B73C", "#2C6DAF", "#40B08F")

tema_slide <- theme_minimal(base_size = 14) +
  theme(panel.grid.minor  = element_blank(),
        plot.background   = element_rect(fill = "transparent", colour = NA),
        panel.background  = element_rect(fill = "transparent", colour = NA),
        legend.position   = "bottom")

# Funcion para exportar cada figura a graphics/ (la carpeta del .tex)
dir.create("graphics", showWarnings = FALSE)
guardar <- function(p, nombre, w = 7, h = 4.2) {
  ggsave(file.path("graphics", paste0(nombre, ".png")), p,
         width = w, height = h, dpi = 300, bg = "transparent")
}

# ---- 2. DATOS -----------------------------------------
data("pottery", package = "HSAUR3")
X     <- pottery[, 1:9]            # solo los 9 oxidos (variables activas)
horno <- factor(pottery$kiln)      # NO se usa para agrupar, solo para validar

dim(X)                             # 45 x 9
summary(X)                         # escalas muy distintas -> hay que estandarizar
round(cor(X), 2)                   # correlaciones altas -> el ACP tiene sentido

# ---- 3. ESTANDARIZAR (z-score) ------------------------
# z_ij = (x_ij - media_j) / s_j   ->  cada variable con media 0 y sd 1
Z <- scale(X)
round(colMeans(Z), 8)              # comprobacion: medias ~ 0
apply(Z, 2, sd)                    # comprobacion: sd = 1

# ---- 4. ACP -------------------------------------------
# scale.unit = TRUE equivale a trabajar con la matriz de correlaciones R
acp <- PCA(X, scale.unit = TRUE, graph = FALSE)
round(acp$eig, 2)                  # valores propios, % varianza y % acumulado

m <- which(acp$eig[, 3] >= 80)[1]  # n. de CP para retener al menos 80% de varianza
m

p_scree <- fviz_eig(acp, addlabels = TRUE, ncp = 9,
                    barfill = col_unal[3], barcolor = col_unal[3],
                    linecolor = col_unal[2], ggtheme = tema_slide) +
  labs(title = NULL, x = "Componente principal", y = "% de varianza")
guardar(p_scree, "01_scree", w = 6, h = 3.6)

scores <- acp$ind$coord[, 1:m]     # puntajes de las m primeras CP

# ---- 5. DISTANCIA + WARD ------------------------------
d  <- dist(scores, method = "euclidean")   # matriz de distancias n x n
hc <- hclust(d, method = "ward.D2")        # OJO: ward.D2, no "ward"

# Etiquetas legibles: horno de origen + numero de pieza
hc$labels <- paste0("H", horno, "-", seq_along(horno))

# Que tan fiel es el arbol a las distancias originales (cercano a 1 = bueno)
cor(d, cophenetic(hc))

# ---- 6. NUMERO DE GRUPOS ------------------------------
# Altura de la fusion que pasa de k a k-1 grupos = rev(hc$height)[k-1]
df_f <- data.frame(k = 2:10, altura = rev(hc$height)[1:9])
df_f$salto <- c(NA, -diff(df_f$altura)) # caida de altura al agregar un grupo

k <- 3   # <-- AJUSTAR despues de mirar el dendrograma y los saltos

p_fus <- ggplot(df_f, aes(k, altura)) +
  geom_line(colour = col_unal[3], linewidth = 1) +
  geom_point(colour = col_unal[3], size = 3) +
  geom_vline(xintercept = k, linetype = 2, colour = col_unal[2]) +
  scale_x_continuous(breaks = 2:10) +
  labs(x = "Numero de grupos k", y = "Altura de fusion (Ward)") +
  tema_slide
guardar(p_fus, "02_fusion", w = 5, h = 3.4)

# ---- 7. DENDROGRAMA -----------------------------------
p_dend <- fviz_dend(hc, k = k, k_colors = col_unal[1:k],
                    rect = TRUE, rect_fill = TRUE, rect_border = "cluster",
                    cex = 0.55, lwd = 0.6, main = "",
                    xlab = "", ylab = "Altura (Ward)",
                    ggtheme = tema_slide)
guardar(p_dend, "03_dendrograma", w = 8, h = 4.6)

# ---- 8. CORTE Y VALIDACION ----------------------------
grupo <- cutree(hc, k = k)         # asigna cada pieza a un grupo
table(grupo)                       # tamano de cada grupo
table(grupo, horno)                # los grupos recuperan los hornos?

# ---- 9. VISUALIZAR LOS GRUPOS EN EL PLANO DEL ACP -----
p_clu <- fviz_cluster(list(data = scores[, 1:2], cluster = grupo),
                      palette = col_unal[1:k], ellipse.type = "convex",
                      repel = TRUE, labelsize = 8, ggtheme = tema_slide) +
  labs(title = NULL,
       x = sprintf("Dim 1 (%.1f%%)", acp$eig[1, 2]),
       y = sprintf("Dim 2 (%.1f%%)", acp$eig[2, 2]))
guardar(p_clu, "04_plano_grupos", w = 6, h = 4.2)

# ---- 10. PERFIL DE CADA GRUPO (media de z por oxido) --
M  <- aggregate(as.data.frame(Z), list(grupo = grupo), mean)
Mm <- as.matrix(M[, -1])
rownames(Mm) <- paste("Grupo", M$grupo)
df_h <- as.data.frame(as.table(Mm))      # columnas: Var1 (grupo), Var2 (oxido), Freq

p_heat <- ggplot(df_h, aes(Var2, Var1, fill = Freq)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(Freq, 2)), size = 3.5) +
  scale_fill_gradient2(low = col_unal[5], mid = "white", high = col_unal[2],
                       midpoint = 0, name = "z") +
  labs(x = NULL, y = NULL) + tema_slide
guardar(p_heat, "05_perfiles", w = 6, h = 3)

# ---- 11. SENSIBILIDAD AL METODO DE ENLACE -------------
metodos <- c("single", "complete", "average", "ward.D2")
sapply(metodos, function(mm) round(cor(d, cophenetic(hclust(d, method = mm))), 3))






















# =====================================================================
# CLUSTERING JERARQUICO (4 criterios de enlace) + ACP  -  dataset pottery
#       y NO han sido ejecutadas aun: corrigelas si algo da error.
# =====================================================================

# ---- 0. PAQUETES --------------------------------------
paquetes <- c("HSAUR3", "FactoMineR", "factoextra", "ggplot2",
              "cluster", "mclust", "patchwork")
nuevos <- paquetes[!paquetes %in% rownames(installed.packages())]
if (length(nuevos) > 0) install.packages(nuevos)
invisible(lapply(paquetes, library, character.only = TRUE))

# ---- 1. ESTILO PARA DIAPOSITIVAS ----------------------
col_unal <- c("#94B43B", "#A61C31", "#35408F", "#F3B73C", "#2C6DAF", "#40B08F")
tema_slide <- theme_minimal(base_size = 14) +
  theme(panel.grid.minor  = element_blank(),
        plot.background   = element_rect(fill = "transparent", colour = NA),
        panel.background  = element_rect(fill = "transparent", colour = NA),
        legend.position   = "bottom")
dir.create("graphics", showWarnings = FALSE)
guardar <- function(p, nombre, w = 7, h = 4.2) {
  ggsave(file.path("graphics", paste0(nombre, ".png")), p,
         width = w, height = h, dpi = 300, bg = "transparent")
}

# ---- 2. DATOS -----------------------------------------
data("pottery", package = "HSAUR3")
X     <- pottery[, 1:9]
horno <- factor(pottery$kiln)          # solo para validar

# ---- 3. ESTANDARIZAR ----------------------------------
Z <- scale(X)

# ---- 4. ACP -------------------------------------------
acp <- PCA(X, scale.unit = TRUE, graph = FALSE)
m   <- which(acp$eig[, 3] >= 80)[1]
p_scree <- fviz_eig(acp, addlabels = TRUE, ncp = 9,
                    barfill = col_unal[3], barcolor = col_unal[3],
                    linecolor = col_unal[2], ggtheme = tema_slide) +
  labs(title = NULL, x = "Componente principal", y = "% de varianza")
guardar(p_scree, "01_scree", w = 6, h = 3.6)
scores <- acp$ind$coord[, 1:m]

# ---- 5. DISTANCIA + WARD ------------------------------
d  <- dist(scores, method = "euclidean")
hc <- hclust(d, method = "ward.D2")
hc$labels <- paste0("H", horno, "-", seq_along(horno))

# ---- 6. NUMERO DE GRUPOS (Ward) -----------------------
df_f <- data.frame(k = 2:10, altura = rev(hc$height)[1:9])
df_f$salto <- c(NA, -diff(df_f$altura))
print(df_f)                            # <- REVISAR antes de fijar k
k <- 3

p_fus <- ggplot(df_f, aes(k, altura)) +
  geom_line(colour = col_unal[3], linewidth = 1) +
  geom_point(colour = col_unal[3], size = 3) +
  geom_vline(xintercept = k, linetype = 2, colour = col_unal[2]) +
  scale_x_continuous(breaks = 2:10) +
  labs(x = "Numero de grupos k", y = "Altura de fusion (Ward)") + tema_slide
guardar(p_fus, "02_fusion", w = 5, h = 3.4)

# ---- 7. DENDROGRAMA (Ward) ----------------------------
p_dend <- fviz_dend(hc, k = k, k_colors = col_unal[1:k],
                    rect = TRUE, rect_fill = TRUE, rect_border = "cluster",
                    cex = 0.55, lwd = 0.6, main = "",
                    xlab = "", ylab = "Altura (Ward)", ggtheme = tema_slide)
guardar(p_dend, "03_dendrograma", w = 8, h = 4.6)

# ---- 8. CORTE Y VALIDACION ----------------------------
grupo <- cutree(hc, k = k)
table(grupo)
table(grupo, horno)

# ---- 9. GRUPOS EN EL PLANO DEL ACP --------------------
p_clu <- fviz_cluster(list(data = scores[, 1:2], cluster = grupo),
                      palette = col_unal[1:k], ellipse.type = "convex",
                      repel = TRUE, labelsize = 8, ggtheme = tema_slide) +
  labs(title = NULL,
       x = sprintf("Dim 1 (%.1f%%)", acp$eig[1, 2]),
       y = sprintf("Dim 2 (%.1f%%)", acp$eig[2, 2]))
guardar(p_clu, "04_plano_grupos", w = 6, h = 4.2)

# ---- 10. PERFIL DE CADA GRUPO -------------------------
M  <- aggregate(as.data.frame(Z), list(grupo = grupo), mean)
Mm <- as.matrix(M[, -1]); rownames(Mm) <- paste("Grupo", M$grupo)
df_h <- as.data.frame(as.table(Mm))
p_heat <- ggplot(df_h, aes(Var2, Var1, fill = Freq)) +
  geom_tile(colour = "white") + geom_text(aes(label = round(Freq, 2)), size = 3.5) +
  scale_fill_gradient2(low = col_unal[5], mid = "white", high = col_unal[2],
                       midpoint = 0, name = "z") +
  labs(x = NULL, y = NULL) + tema_slide
guardar(p_heat, "05_perfiles", w = 6, h = 3)

# ---- 11. COFENETICA (tu salida previa) ----------------
metodos <- c(Simple = "single", Completo = "complete",
             Promedio = "average", Ward = "ward.D2")
hcs  <- lapply(metodos, function(mm) hclust(d, method = mm))
coph <- sapply(hcs, function(h) cor(d, cophenetic(h)))
round(coph, 3)

# =====================================================================
# NUEVO 12. COMPARAR LOS 4 CRITERIOS: DENDROGRAMAS LADO A LADO
# =====================================================================
dends <- lapply(names(hcs), function(nm) {
  h <- hcs[[nm]]; h$labels <- hc$labels
  fviz_dend(h, k = k, k_colors = col_unal[1:k], cex = 0.3, lwd = 0.5,
            main = nm, xlab = "", ylab = "Altura", ggtheme = theme_minimal())
})
p_cuatro <- wrap_plots(dends, ncol = 2)
guardar(p_cuatro, "06_cuatro_enlaces", w = 10, h = 6.5)

# =====================================================================
# NUEVO 13. CRITERIOS PARA ELEGIR EL ENLACE Y k
#   (a) correlacion cofenetica   (b) coeficiente aglomerativo (agnes)
#   (c) silueta media por k      (d) acuerdo con la validacion externa (ARI)
# =====================================================================
ag_met <- c(Simple = "single", Completo = "complete",
            Promedio = "average", Ward = "ward")   # nombre en agnes()
ac <- sapply(ag_met, function(mm) agnes(d, method = mm)$ac)

sil <- sapply(hcs, function(h)
  sapply(2:8, function(kk) mean(silhouette(cutree(h, kk), d)[, "sil_width"])))
rownames(sil) <- 2:8

ari_horno <- sapply(hcs, function(h) adjustedRandIndex(cutree(h, k), horno))

resumen <- data.frame(Cofenetica = round(coph, 3),
                      Aglomerativo = round(ac, 3),
                      Silueta_k = round(sil[as.character(k), ], 3),
                      ARI_horno = round(ari_horno, 3))
print(resumen)                         # <- tabla para la diapositiva

df_s <- data.frame(k = rep(2:8, times = 4),
                   metodo = rep(names(hcs), each = 7),
                   silueta = as.vector(sil))
p_sil <- ggplot(df_s, aes(k, silueta, colour = metodo)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = col_unal[c(3, 2, 1, 4)], name = NULL) +
  scale_x_continuous(breaks = 2:8) +
  labs(x = "Numero de grupos k", y = "Silueta media") + tema_slide
guardar(p_sil, "07_silueta", w = 6, h = 3.6)

# Silueta individual de Ward con k elegido
p_sil_w <- fviz_silhouette(silhouette(grupo, d), palette = col_unal[1:k],
                           ggtheme = tema_slide, print.summary = FALSE) +
  labs(title = NULL)
guardar(p_sil_w, "08_silueta_ward", w = 6, h = 3.6)

# =====================================================================
# NUEVO 14. CONSISTENCIA: repetir codigo, remuestreo y comparacion con k-means
# =====================================================================
# (a) Repetir el codigo: hclust es DETERMINISTA (no usa azar)
hc2 <- hclust(dist(scores), method = "ward.D2")
identical(unname(cutree(hc, k)), unname(cutree(hc2, k)))          # TRUE

# k-means SI depende de la semilla de inicio (1 inicio por corrida)
km_ref  <- kmeans(scores, centers = k, nstart = 1)$cluster
ari_km  <- sapply(1:30, function(s) { set.seed(s)
  adjustedRandIndex(km_ref, kmeans(scores, centers = k, nstart = 1)$cluster) })

# (b) Otra muestra de la poblacion: bootstrap de las piezas
set.seed(123)
B <- 200
ari_boot <- replicate(B, {
  idx <- sample(nrow(scores), replace = TRUE)
  g_b <- cutree(hclust(dist(scores[idx, ]), method = "ward.D2"), k)
  adjustedRandIndex(grupo[idx], g_b)
})
summary(ari_boot)                       # cerca de 1 = grupos estables

df_est <- rbind(data.frame(prueba = "Ward bootstrap", ARI = ari_boot),
                data.frame(prueba = "K-means semillas", ARI = ari_km))
p_est <- ggplot(df_est, aes(prueba, ARI, fill = prueba)) +
  geom_boxplot() + scale_fill_manual(values = col_unal[c(1, 4)], guide = "none") +
  coord_cartesian(ylim = c(0, 1)) + labs(x = NULL, y = "ARI vs. referencia") + tema_slide
guardar(p_est, "09_estabilidad", w = 5, h = 3.4)

# =====================================================================
# NUEVO 15. ATIPICOS: agregar una pieza extrema y ver que pasa
# =====================================================================
set.seed(1)
out    <- matrix(c(6, 6, 6)[1:m], nrow = 1)        # punto muy lejano en las m CP
sc_out <- rbind(scores, out); rownames(sc_out)[nrow(sc_out)] <- "OUT"
d_out  <- dist(sc_out)

res_out <- sapply(metodos, function(mm) {
  g <- cutree(hclust(d_out, method = mm), k)
  c(ARI_originales = round(adjustedRandIndex(grupo, g[1:45]), 3),
    tam_min        = min(table(g)))                # tam_min = 1 -> el outlier es "su grupo"
})
print(res_out)

# Alternativa resistente: k-medoides (PAM) usa medoides en vez de medias
pam_out <- pam(d_out, k = k)
adjustedRandIndex(grupo, pam_out$clustering[1:45])

# Deteccion previa: distancia de Mahalanobis al centro (libro, cap. 2)
d2 <- mahalanobis(scores, colMeans(scores), cov(scores))
which(d2 > qchisq(0.975, df = m))       # candidatas a atipicas

# Las 5 observaciones más alejadas del centro
sort(d2, decreasing = TRUE)[1:5]

# ¿Qué porcentaje supera el umbral?
mean(d2 > qchisq(0.975, df = m)) * 100

# =====================================================================
# NUEVO 16. ABREBOCAS: alternativa basada en modelos (mixturas gaussianas)
# =====================================================================
mc <- Mclust(scores)
summary(mc)                             # elige k y forma por BIC
adjustedRandIndex(mc$classification, horno)

# =====================================================================
# BONUS SHINY (exploracion rapida)
# =====================================================================
# library(Factoshiny)
# res.pca  <- PCAshiny(X)
# res.hcpc <- HCPCshiny(res.pca)