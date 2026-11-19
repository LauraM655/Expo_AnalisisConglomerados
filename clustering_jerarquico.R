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
