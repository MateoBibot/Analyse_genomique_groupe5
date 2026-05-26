library(ape)

dna <- read.dna("data/ASFV_genomic_analyses_simulated_dataset_1-2.fas", format = "fasta")

# Raccourcir les noms à 6 caractères max (obligation de Network)
noms_originaux <- rownames(dna)
noms_courts    <- make.unique(substr(noms_originaux, 1, 6), sep = "")
rownames(dna)  <- noms_courts

output_dir <- "./output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
output_path <- file.path(output_dir, "ASFV_pour_Network_v2.nex")

# Exporter en NEXUS séquentiel
write.nexus.data(dna, 
                 file = output_path,  
                 format = "dna",
                 interleaved = FALSE)

cat("Fichier créé :", output_path, "\n")

# Garde la correspondance noms courts <-> noms originaux pour plus tard
correspondance <- data.frame(nom_court    = noms_courts,
                             nom_original = noms_originaux)
write.csv(correspondance, "correspondance_noms.csv", row.names = FALSE)




# MODÉLISATION SPATIO-TEMPORELLE DE LA PROPAGATION DE LA FIÈVRE PORCINE AFRICAINE (ASFV)
setwd("C:/Users/mbibo/OneDrive - Université Libre de Bruxelles/BING-F432/nouvelle version") #à modifier avec votre propre chemin d'accès
# --- Chargement des packages nécessaires ---

if (!require(fields)) install.packages("fields")
if (!require(raster)) install.packages("raster")
if (!require(sf)) install.packages("sf")
if (!require(sp)) install.packages("sp")
if (!require(ks)) install.packages("ks")
if (!require(RColorBrewer)) install.packages("RColorBrewer")
if (!require(lubridate)) install.packages("lubridate")
if(!require(EpiEstim))install.packages("EpiEstim")
library(fields); library(raster); library(sf); library(sp); library(ks); library(RColorBrewer); library(lubridate);library(EpiEstim);


#estimation du nombre de reproduction effectif (R(t))
data0 = read.csv("ASFV_genomic_analyses_simulated_dataset_1-2.csv", header = TRUE)
head(data0)

# Calcul des jours écoulés et de la durée totale de l'étude
days = interval(min(ymd(data0[, "collection_date"])), ymd(data0[, "collection_date"])) %/% days(1) + 1
total_number_of_days = interval(min(ymd(data0[, "collection_date"])), max(ymd(data0[, "collection_date"]))) %/% days(1) - 1

daily_cases = rep(NA, total_number_of_days)
for(i in 1:length(daily_cases)) {
  daily_cases[i] = sum(days == i)
}

n = 1000 # Nombre d'itérations

# --- MODIFICATION IMPORTANTE SUITE AU TP ---
mean_range = c(13, 18) # Plage de la moyenne de l'intervalle sériel (SI) spécifique à ASFV
#L'intervalle sériel est le temps entre l'infection d'un hôte et l'infection de ses contacts.
sd_range = c(5, 7)     # Plage de l'écart-type du SI spécifique à ASFV
# -------------------------------------------

# On calcule le R(t) sur des blocs de 7 jours pour lisser les variations quotidiennes
t_start = seq(2, length(daily_cases) - 6) # Fenêtre glissante de 7 jours
t_end = seq(8, length(daily_cases))
all_Rt = matrix(NA, nrow = length(t_start), ncol = n) # Matrice pour stocker les estimations

# Estimation de Rt sur n itérations avec tirage uniforme
for(i in 1:n) {
  mean_si_i = runif(1, mean_range[1], mean_range[2])
  sd_si_i = runif(1, sd_range[1], sd_range[2])
  
  res_i = estimate_R(incid = daily_cases, method = "parametric_si", config = make_config(list(
    mean_si = mean_si_i, std_si = sd_si_i, t_start = t_start, t_end = t_end
  )))
  all_Rt[, i] = res_i$R$`Mean(R)`
}

# Calcul de la médiane, des dates correspondantes et des intervalles de confiance
R_median = apply(all_Rt, 1, median, na.rm = TRUE)
Rt_days = (t_start + t_end) / 2
R_dates = decimal_date(min(ymd(data0[, "collection_date"])) + Rt_days)

R_lower = apply(all_Rt, 1, quantile, probs = 0.025, na.rm = TRUE)
R_upper = apply(all_Rt, 1, quantile, probs = 0.975, na.rm = TRUE)

# --- Configuration et tracé du graphique ---

# Tracé principal : ajout de xaxs="i", yaxs="i" et modification de xlim
par(oma=c(0,0,0,0),mar=c(2.0,3.5,1,2),lwd=0.3,bty="o",col="gray30",col.axis="gray30",
    fg="gray30")
plot(R_dates, R_median, lwd = 0.7, type = "l", cex.axis = 0.8, cex.lab = 0.8, col = "gray30", 
     axes = FALSE, xlab = NA, ylab = NA, 
     ylim = c(0, 18), 
     xlim = c(decimal_date(ymd("2024-02-01")), 2025.0), # Début en février
     xaxs = "i", yaxs = "i") 

# Bande de confiance (polygone gris)
xx_l = c(R_dates, rev(R_dates))
yy_l = c(R_lower, rev(R_upper))
polygon(xx_l, yy_l, col = rgb(187/255, 187/255, 187/255, 0.25), border = 0)

# Ligne médiane et seuil de 1
lines(R_dates, R_median, lwd = 0.3, type = "l", cex.axis = 0.8, cex.lab = 0.8, col = "gray30")
abline(h = 1, lty = 2, lwd = 0.3, col = "gray30")

# Définition des dates pour l'axe des abscisses 
dates = c("2024-02-01", "2024-03-01", "2024-04-01", "2024-05-01", "2024-06-01", 
          "2024-07-01", "2024-08-01", "2024-09-01", "2024-10-01", "2024-11-01", "2024-12-01", "2025-01-01")
ats = decimal_date(ymd(dates))
labels = format(ymd(dates), "%d/%m/%Y") 


# Création des axes
axis(side = 1, lwd = 0.5, cex.axis = 0.7, mgp = c(0, 0.5, 0), lwd.tick = 0.5, col.lab = "gray30", 
     col = "gray30", tck = -0.03, las = 1, at = ats, label = labels)
axis(side = 2, lwd = 0.5, cex.axis = 0.7, mgp = c(0, 0.7, 0), lwd.tick = 0.5, col.lab = "gray30", 
     col = "gray30", tck = -0.03, las = 1, padj = 0.4, at = seq(0, 18, by = 2)) # Ajustement des ticks Y
# Légende de l'axe Y
mtext("Taux de reproduction effectif", side = 2, col = "gray30", cex = 0.9, line = 1.7, las = 3)



# --- Chargement des couches géographiques ---

borders = shapefile("data/Country_borders_in_study_area.shp")
motorways = shapefile("data/Motorways_lines_in_study_area.shp")
land_cover_variables = raster("data/Raster_Wallonian_forest_areas.asc")

crs_lb72 = CRS("+init=epsg:31370")

# On assigne le CRS au raster
crs(land_cover_variables) = crs_lb72

# On reprojette les shapefiles vers ce même CRS
borders_r   = spTransform(borders,   crs_lb72)
motorways_r = spTransform(motorways, crs_lb72)



# Chaque ligne correspond à un cas détecté d'ASFV, avec ses coordonnées et sa date de collecte.
data1 = read.csv("data/ASFV_genomic_analyses_simulated_dataset_1-2.csv", header=TRUE,)
head(data1)


# On trie du plus ancien au plus récent avant tout calcul.
data1 = data1[order(data1$collection_date), ]

# Calcul du nombre de jours écoulés depuis le premier cas détecté.
# La fonction interval() de lubridate mesure la durée entre deux dates,
# et /days(1) convertit cet intervalle en nombre de jours (valeur numérique).
data1$days_collection <- interval(data1$collection_date[1], data1$collection_date)/days(1)

head(data1)



# On utilise le raster de couverture terrestre comme gabarit:ses dimensions  définissent la grille sur laquelle l'estimation par noyau sera calculée.
rast = land_cover_variables
gridSize = c(rast@ncols, rast@nrows) # nombre de colonnes et lignes de la grille
xyMin = c(rast@extent@xmin, rast@extent@ymin); xyMax = c(rast@extent@xmax, rast@extent@ymax)

# La bande passante H est super importante pour l'estimation par noyau.Elle contrôle le lissage de la densité estimée. 
H = ks::Hpi(data1[,c("longitude","latitude")])

# data2 stockera les cas filtrés (front d'onde), buffer accumule tous les cas  pour recalculer le KDE à chaque itération.
data2 = data1[which(data1[,"days_collection"]==0),]
buffer = data1[which(data1[,"days_collection"]==0),]

# Création du dossier de sortie si inexistant.
if (!dir.exists("KDE95_contour_polygons")) {
  dir.create("KDE95_contour_polygons")
}

# L'objectif est de ne conserver que les cas situés au "front d'onde" de l'épidémie, c'est-à-dire les cas qui étendent le périmètre infecté.
# Les cas apparaissant à l'intérieur de la zone déjà contaminée (contour 95% du KDE) sont considérés comme des réinfections locales et sont écartés.
for (i in 1:max(data1[,"days_collection"])) {
  cat("\tDay",i,"\n",sep=" ")
  indices1 = which(data1[,"days_collection"]==i)
  if (length(indices1) > 0) { # on ne recalcule le KDE que si de nouveaux cas apparaissent ce jour-là
    if (dim(data2)[1] >= 3) { # le KDE nécessite au minimum 3 points pour être calculé
      
      
      kde = ks::kde(buffer[,c("longitude","latitude")], compute.cont=T,
                    H=H, gridsize=gridSize, xmin=xyMin, xmax=xyMax)
      
      # Conversion du KDE en raster, puis extraction du contour à 95%.
      # Le contour à 95% délimite la région qui contient 95% de la probabilité de présence estimée
      r = raster(kde)
      contour = rasterToContour(r, levels=kde$cont["5%"])
      
      # On identifie le seuil de valeur raster correspondant au contour 95%,
      # puis on met à NA toutes les cellules en dessous (hors de la zone infectée).
      threshold = contourLevels(kde, 0.05)
      r[r[]<threshold] = NA
      r[!is.na(r[])] = i # les cellules dans la zone infectée reçoivent le numéro du jour correspondant
      crs(r) = crs(rast)
      
      # Sauvegarde du raster et du shapefile du contour 95% pour ce jour.
      file_name = paste0("KDE95_contour_polygons/KDE_95_contour_day_",i-1,".tif")
      writeRaster(r, file_name, overwrite=T)
      file_name = paste0("KDE95_contour_polygons/KDE_95_contour_day_",i-1,".shp")
      st_write(st_as_sf(contour), file_name, update=T)
      
      # Pour chaque nouveau cas du jour i, on vérifie s'il se trouve à l'intérieur du contour 95% déjà établi. Si oui, il est ignoré (cas interne).
      # Si non, il est conservé car il repousse le front d'onde (cas externe).
      indices2 = c()
      for (j in 1:length(indices1)) {
        point_in_polygon = FALSE
        pt.x = data1[indices1[j],"longitude"]
        pt.y = data1[indices1[j],"latitude"]
        for (k in 1:length(contour@lines[[1]]@Lines)) {
          pol.x = contour@lines[[1]]@Lines[[k]]@coords[,1]
          pol.y = contour@lines[[1]]@Lines[[k]]@coords[,2]
          # point.in.polygon() retourne 0 si le point est hors du polygone, autre valeur sinon
          if (point.in.polygon(pt.x, pt.y, pol.x, pol.y) != 0) {
            point_in_polygon = TRUE
          }
        }
        if (point_in_polygon == FALSE) {
          indices2 = c(indices2, indices1[j]) # cas conservé car il étend le front d'onde
        }
      }
      if (length(indices2) > 0) {
        data2 = rbind(data2, data1[indices2,])
      }
    } else {
      
      data2 = rbind(data2, data1[indices1,])
    }
    buffer = rbind(buffer, data1[indices1,]) # on accumule tous les cas (filtrés ou non) dans le buffer
  }
}

# Suppression des doublons 
data2 = unique(data2)
write.csv(data2, "ASFV_filtered_cases.csv", quote=F, row.names=F)



# On fusionne successivement chaque raster KDE journalier sauvegardé.
buffer = land_cover_variables; buffer[] = NA


for (i in 1:max(data1[,"days_collection"])) {
  if (file.exists(paste0("KDE95_contour_polygons/KDE_95_contour_day_",i,".tif"))) {
    kde = raster(paste0("KDE95_contour_polygons/KDE_95_contour_day_",i,".tif"))
    buffer = merge(buffer, kde)
  }
}
# Sauvegarde du raster cumulatif final et visualisation rapide.
writeRaster(buffer, "KDE95_contour_raster.tif", overwrite=T)
plot(buffer)


# On estime les cas d'infection de façon continue (nuage de densité) plutôt que de travailler avec des points discrets. A ~350 jours, les cas
# les plus récents sont à la périphérie du front d'onde.
# L'objectif est de créer un masque géographique qui délimite la zone réellement touchée par l'épidémie. 

data2 = read.csv("ASFV_filtered_cases.csv", head=T)

# On crée un raster gabarit avec des zéros partout où la donnée n'est pas NA.
template = land_cover_variables; template[!is.na(template[])] = 0

H = Hpi(data2[,c("longitude","latitude")])
kde = kde(data2[,c("longitude","latitude")],
          H=H, compute.cont=T, gridsize=c(1000,1000))
rast1 = raster(kde); contour = rasterToContour(rast1, levels=kde$cont["5%"])
coords = data2[,c("longitude","latitude")]
threshold = min(raster::extract(rast1, coords))

p = Polygon(contour@lines[[1]]@Lines[[1]]@coords)
ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
rast2 = mask(rast1, sps); mask = raster::resample(rast2, template)

# On veut estimer le premier temps d'invasion en tout point de l'espace à partir des cas filtrés, en produisant une surface lisse et continue.
coords = data2[,c("longitude","latitude")]
tps_model = Tps(x=coords, Y=data2[,"days_collection"]) # entraînement du modèle TPS
tps = interpolate(template, tps_model) # interpolation sur toute la grille du gabarit
tps_mask = mask(tps, mask) # application du masque : on ne garde que la zone envahie

plot(tps_mask)


# La friction mesure la vitesse locale de changement du temps d'invasion entre une cellule et ses 8 voisines. Une forte variation temporelle entre cellules
# adjacentes indique une propagation lente (haute friction).
raster_resolution = mean(c(res(tps_mask)[1],res(tps_mask)[2]))
f = matrix(1/raster_resolution, nrow=3, ncol=3)
f[c(1,3,7,9)] = 1/(sqrt(2)*raster_resolution); f[5] = 0 # poids nul pour la cellule centrale
fun = function(x, ...) {
  sum(abs(x-x[5])*f)/8 # différence moyenne pondérée entre la cellule centrale et ses 8 voisines
}
friction = focal(tps_mask, w=matrix(1,nrow=3,ncol=3), fun=fun, pad=T, padValue=NA, na.rm=F)


myAvSize = 11
friction_sd10 = focal(friction,
                      w=matrix(1/(myAvSize^2), nrow=myAvSize, ncol=myAvSize),
                      pad=T, padValue=NA, na.rm=F)


# La vitesse du front d'onde est l'inverse de la friction (une friction faible = propagation rapide). On divise par 1000 pour convertir les mètres en kilomètres,
# puis on multiplie par 7 pour passer du jour à la semaine.
wavefrontVelocity_sd10 = ((1/(friction_sd10))/1000)*7
print(mean(wavefrontVelocity_sd10[], na.rm=T)) 

#--------------------------------------------------------------------
crs_ref = crs(land_cover_variables)   #on garde le CRS du raster

# Reprojection des couches vectorielles
borders_r   = spTransform(borders,   crs_ref)
motorways_r = spTransform(motorways, crs_ref)
cat("Extent raster   :", as.character(extent(land_cover_variables)), "\n")
cat("Extent motorways:", as.character(extent(motorways_r)), "\n")
cat("Extent borders  :", as.character(extent(borders_r)), "\n")


# On passe maintenant à la visualisation des données



library(raster); library(sp); library(sf)
library(RColorBrewer); library(fields)


data1$x_proj = data1$longitude
data1$y_proj = data1$latitude
data2$x_proj = data2$longitude
data2$y_proj = data2$latitude

crs_raster = CRS("EPSG:3035")
crs(land_cover_variables) = crs_raster
borders_r   = spTransform(borders,   crs_raster)
motorways_r = spTransform(motorways, crs_raster)




cols_palette = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(161)[21:121])

cols_lc = c(
  rgb(0, 1, 0, 0),
  rgb(144/255, 238/255, 144/255, 0.4),
  rgb(0, 1, 0, 0),
  rgb(190/255, 190/255, 190/255, 0.6)
)

col_motorway = rgb(220/255, 50/255, 30/255, 1)
col_border   = "black"




add_colorbar = function(pal, zlim, at_vals, labels=NULL,
                        width=0.018, shrink=0.70, gap=0.012) {
  usr    = par("usr")
  plot_w = usr[2] - usr[1]
  plot_h = usr[4] - usr[3]
  x0 = usr[2] + gap * plot_w
  x1 = usr[2] + (gap + width) * plot_w
  margin_v = (1 - shrink) / 2
  y0 = usr[3] + margin_v * plot_h
  y1 = usr[4] - margin_v * plot_h
  n   = length(pal)
  y_s = seq(y0, y1, length.out = n + 1)
  for (k in seq_len(n)) rect(x0, y_s[k], x1, y_s[k+1], col=pal[k], border=NA, xpd=TRUE)
  rect(x0, y0, x1, y1, col=NA, border="gray40", lwd=0.4, xpd=TRUE)
  if (is.null(labels)) labels = as.character(at_vals)
  at_y = y0 + (at_vals - zlim[1]) / diff(zlim) * (y1 - y0)
  for (k in seq_along(at_vals)) {
    segments(x1, at_y[k], x1 + 0.005*plot_w, at_y[k], xpd=TRUE, lwd=0.4, col="gray40")
    text(x1 + 0.008*plot_w, at_y[k], labels=labels[k],
         xpd=TRUE, cex=0.60, adj=c(0, 0.5), col="gray20")
  }
}



draw_vectors = function(lwd_mway=1.8, lwd_border=0.9) {
  # Pour réafficher les frontières et routes, enlevez simplement les '#' devant les 4 lignes 'plot' ci-dessous.
  
   plot(borders_r,   add=TRUE, lwd=lwd_border * 3, col="white",       lty=1)
   plot(borders_r,   add=TRUE, lwd=lwd_border,      col=col_border,    lty=2)
   plot(motorways_r, add=TRUE, lwd=lwd_mway * 2.5, col="white",       lty=1)
   plot(motorways_r, add=TRUE, lwd=lwd_mway,        col=col_motorway,  lty=1)
}


# Mise en page 


par(mfrow = c(2, 2),
    mar   = c(0.8, 0.8, 1.8, 5.5),
    oma   = c(0.5, 0, 2.5, 0),
    mgp   = c(0, 0.4, 0),
    lwd   = 0.3,
    bty   = "o")


# PANNEAU 1 — Tous les cas

days1     = data1[, "days_collection"]
cols_pts1 = cols_palette[(((days1 - min(days1)) / (max(days1) - min(days1))) * 100) + 1]

# 1. Fond raster
plot(land_cover_variables, col=cols_lc, box=FALSE, axes=FALSE, legend=FALSE)

# 2. Points (anciens par-dessus)
for (i in nrow(data1):1) {
  points(data1[i, "x_proj"], data1[i, "y_proj"], pch=16, cex=0.65, col=cols_pts1[i])
  points(data1[i, "x_proj"], data1[i, "y_proj"], pch=1,  cex=0.65, col="gray30", lwd=0.3)
}

# 3. Vecteurs par dessus les points 
draw_vectors()

box(lwd=0.4, col="gray30")
mtext("1. Premiers temps d'invasion", side=3, line=0.5,  adj=0.05, cex=0.75, font=1)
mtext("(tous, en jours)",             side=3, line=-0.3, adj=0.05, cex=0.65, col="gray30")

at_d = seq(0, max(days1), by=30)
add_colorbar(cols_palette, c(0, max(days1)), at_d)


# PANNEAU 2 — Cas filtrés + légende

days2     = data2[, "days_collection"]
cols_pts2 = cols_palette[(((days2 - min(days2)) / (max(days2) - min(days2))) * 100) + 1]

plot(land_cover_variables, col=cols_lc, box=FALSE, axes=FALSE, legend=FALSE)

for (i in nrow(data2):1) {
  points(data2[i, "x_proj"], data2[i, "y_proj"], pch=16, cex=0.65, col=cols_pts2[i])
  points(data2[i, "x_proj"], data2[i, "y_proj"], pch=1,  cex=0.65, col="gray30", lwd=0.3)
}

draw_vectors()

box(lwd=0.4, col="gray30")
mtext("2. Premiers temps d'invasion", side=3, line=0.5,  adj=0.05, cex=0.75, font=1)
mtext("(filtrés, en jours)",          side=3, line=-0.3, adj=0.05, cex=0.65, col="gray30")

 legend("bottomleft",
        legend   = c("Autoroutes", "Frontière int."),
        col      = c(col_motorway, col_border),
        lty      = c(1, 2),
        lwd      = c(1.8, 0.9),
        cex      = 0.62,
        bty      = "n",
        inset    = c(0.01, 0.02),
        text.col = "gray20",
        seg.len  = 1.5)

at_d2 = seq(0, max(days2), by=30)
add_colorbar(cols_palette, c(0, max(days2)), at_d2)



# PANNEAU 3 — Temps interpolés


plot(tps_mask, main="", box=FALSE, axes=FALSE, legend=FALSE,
     col=cols_palette, colNA="white")
plot(land_cover_variables, add=TRUE, col=cols_lc, box=FALSE, axes=FALSE, legend=FALSE)

# Croix des cas filtrés
points(data2[, "x_proj"], data2[, "y_proj"], pch=3, cex=0.45, lwd=0.4, col="gray30")

# Vecteurs par-dessus
draw_vectors()

box(lwd=0.4, col="gray30")
mtext("3. Premiers temps d'invasion", side=3, line=0.5,  adj=0.05, cex=0.75, font=1)
mtext("(interpolés, en jours)",       side=3, line=-0.3, adj=0.05, cex=0.65, col="gray30")

max_tps = max(tps_mask[], na.rm=TRUE)
at_tps  = seq(0, floor(max_tps / 30) * 30, by=30)
add_colorbar(cols_palette, c(0, max_tps), at_tps)



# PANNEAU 4 — Vitesse du front d'onde


zlim_vel    = c(0.1, 1.0)
# Assurez-vous que la variable `wavefrontVelocity_sd10_truncated` est définie plus tôt dans votre code complet. 
# Si ce n'est pas le cas, on utilise ici `wavefrontVelocity_sd10`.
vel_display = wavefrontVelocity_sd10
vel_display[vel_display[] < zlim_vel[1]] = NA

plot(vel_display, main="", box=FALSE, axes=FALSE, legend=FALSE,
     col=cols_palette, colNA="white", zlim=zlim_vel)
plot(land_cover_variables, add=TRUE, col=cols_lc, box=FALSE, axes=FALSE, legend=FALSE)

points(data2[, "x_proj"], data2[, "y_proj"], pch=3, cex=0.45, lwd=0.4, col="gray30")

draw_vectors()

box(lwd=0.4, col="gray30")
mtext("4. Vitesse du front d'onde", side=3, line=0.5,  adj=0.05, cex=0.75, font=1)
mtext("(km/semaine)",               side=3, line=-0.3, adj=0.05, cex=0.65, col="gray30")

at_vel  = seq(0.1, 1.0, by=0.1)
lab_vel = c(as.character(seq(0.1, 0.9, by=0.1)), ">=1")
add_colorbar(cols_palette, zlim_vel, at_vel, lab_vel)



# Titre global

mtext("ASFV — Progression du front d'onde — Sud de la Belgique",
      outer=TRUE, side=3, line=1.0, cex=0.90, font=2, col="gray20")






# --- Création de l'animation visuelle des cas ---

# Chargement des librairies spécifiques à l'animation
required_pkgs = c("gifski", "png", "RColorBrewer", "raster", "sp", "viridisLite")
for (pkg in required_pkgs) {
  if (!require(pkg, character.only=TRUE)) install.packages(pkg)
  library(pkg, character.only=TRUE)
}

# Paramétrage des projections
crs_raster = CRS("EPSG:3035")
crs(land_cover_variables) = crs_raster
borders_r   = spTransform(borders,   crs_raster)
motorways_r = spTransform(motorways, crs_raster)

data1$x_proj = data1$longitude
data1$y_proj = data1$latitude

# Préparation du dossier d'images temporaires
frames_dir = "animation_frames_v2"
if (!dir.exists(frames_dir)) dir.create(frames_dir)

# Paramètres de temps et de couleurs pour l'animation
max_day  = max(data1[, "days_collection"])
days_seq = seq(5, max_day, by=3)

cols_anim = viridis(101, option="plasma")

cols_lc_dark = c(
  rgb(0, 0, 0, 0),
  rgb(30/255, 80/255, 30/255, 0.55),
  rgb(0, 0, 0, 0),
  rgb(35/255, 35/255, 40/255, 1)
)

bg_col    = "#12121a"
col_kde   = rgb(1, 0.35, 0.05, 0.30)
col_halo  = rgb(1, 0.95, 0.3, 0.85)
col_mway   = rgb(1,    0.35, 0.15, 1.0)   
col_border = rgb(0.85, 0.85, 1.0,  0.85)  

# Fonction de dessin de la barre de progression
draw_progress_bar = function(xfrac, label_day, max_day,
                             col_fill="#f0c040", col_bg="#2a2a3a", col_text="white") {
  usr = par("usr")
  w = usr[2] - usr[1]; h = usr[4] - usr[3]
  x0 = usr[1] + 0.04*w;  x1 = usr[2] - 0.04*w
  y0 = usr[3] + 0.012*h; y1 = usr[3] + 0.032*h
  rect(x0, y0, x1, y1, col=col_bg, border=NA)
  rect(x0, y0, x0 + xfrac*(x1-x0), y1, col=col_fill, border=NA)
  rect(x0, y0, x1, y1, col=NA, border=rgb(1,1,1,0.3), lwd=0.5)
  text((x0+x1)/2, (y0+y1)/2,
       labels=sprintf("Jour %d  /  %d   |   Semaine épidémique %d",
                      label_day, max_day, ceiling(label_day/7)),
       col=col_text, cex=0.52, font=2)
}

# Fonction d'attribution des couleurs selon les jours
day_to_col = function(d, max_d, palette) {
  idx = round((d / max_d) * 100) + 1
  idx = pmax(1, pmin(length(palette), idx))
  palette[idx]
}

# Fonction de dessin des frontières et routes pour le GIF
draw_vectors_gif = function() {
  # Pour réafficher les frontières et routes sur le GIF, enlevez simplement les '#' devant les 4 lignes 'plot' ci-dessous.
  
   plot(borders_r,   add=TRUE, lwd=2.5, col=rgb(0,0,0,0.8),  lty=1)
   plot(borders_r,   add=TRUE, lwd=0.9, col=col_border,       lty=2)
   plot(motorways_r, add=TRUE, lwd=4.0, col=rgb(0,0,0,0.7),  lty=1)
   plot(motorways_r, add=TRUE, lwd=1.8, col=col_mway,         lty=1)
}

# Génération des images individuelles (frames)
cat("Génération de", length(days_seq), "frames...\n")

for (i in seq_along(days_seq)) {
  d       = days_seq[i]
  xfrac   = d / max_day
  sub     = data1[data1$days_collection <= d, ]
  days_sub = sub$days_collection
  n_cases = nrow(sub)
  cols_pts = day_to_col(days_sub, max_day, cols_anim)
  is_recent = (d - days_sub) <= 15
  
  kde_file = paste0("KDE95_contour_polygons/KDE_95_contour_day_", d, ".tif")
  kde_shp  = paste0("KDE95_contour_polygons/KDE_95_contour_day_", d, ".shp")
  
  frame_path = sprintf("%s/frame_%04d.png", frames_dir, i)
  png(frame_path, width=1200, height=900, res=150, bg=bg_col)
  par(mar=c(2.5, 0.5, 3.2, 5.5), bg=bg_col, col.main="white", col.lab="white")
  
  # Fond raster
  plot(land_cover_variables, col=cols_lc_dark, box=FALSE, axes=FALSE, legend=FALSE)
  
  # Zone KDE 95%
  if (file.exists(kde_file)) {
    kde_r = raster(kde_file)
    kde_r[!is.na(kde_r[])] = 1
    plot(kde_r, add=TRUE, col=col_kde, legend=FALSE)
    if (file.exists(kde_shp)) {
      kde_border = shapefile(kde_shp)
      plot(kde_border, add=TRUE, lwd=0.8, border=rgb(1,0.55,0.1,0.7), col=NA)
    }
  }
  
  # Points anciens
  if (any(!is_recent)) {
    points(sub[!is_recent, "x_proj"], sub[!is_recent, "y_proj"],
           pch=16, cex=0.55, col=cols_pts[!is_recent])
    points(sub[!is_recent, "x_proj"], sub[!is_recent, "y_proj"],
           pch=1,  cex=0.55, col=rgb(1,1,1,0.15), lwd=0.3)
  }
  
  # Points récents (front actif)
  if (any(is_recent)) {
    points(sub[is_recent, "x_proj"], sub[is_recent, "y_proj"],
           pch=16, cex=2.0, col=rgb(1,0.9,0.2,0.18))
    points(sub[is_recent, "x_proj"], sub[is_recent, "y_proj"],
           pch=16, cex=1.2, col=rgb(1,0.9,0.2,0.35))
    points(sub[is_recent, "x_proj"], sub[is_recent, "y_proj"],
           pch=16, cex=0.75, col=cols_pts[is_recent])
    points(sub[is_recent, "x_proj"], sub[is_recent, "y_proj"],
           pch=1,  cex=0.75, col=rgb(1,1,1,0.7), lwd=0.5)
  }
  
  # Application des tracés vectoriels
  draw_vectors_gif()
  
  # Ajout des titres sur l'image
  usr = par("usr")
  rect(usr[1], usr[4] - 0.055*(usr[4]-usr[3]),
       usr[2], usr[4], col=rgb(0,0,0,0.55), border=NA, xpd=TRUE)
  mtext("ASFV — Propagation du front d'onde épidémique",
        side=3, line=1.4, cex=0.95, font=2, col="white", adj=0.02)
  mtext(sprintf("n = %d cas cumulés  |  Wallonie, Belgique", n_cases),
        side=3, line=0.3, cex=0.68, font=1, col=rgb(1,1,1,0.65), adj=0.02)
  
  # Affichage des barres d'informations
  draw_progress_bar(xfrac, label_day=d, max_day=max_day)
  
  n_leg    = 50
  leg_cols = day_to_col(seq(0, max_day, length.out=n_leg), max_day, cols_anim)
  x_leg0   = usr[2] + 0.01*(usr[2]-usr[1])
  x_leg1   = usr[2] + 0.025*(usr[2]-usr[1])
  y_vals   = seq(usr[3] + 0.1*(usr[4]-usr[3]),
                 usr[4] - 0.1*(usr[4]-usr[3]),
                 length.out=n_leg+1)
  for (k in 1:n_leg)
    rect(x_leg0, y_vals[k], x_leg1, y_vals[k+1], col=leg_cols[k], border=NA, xpd=TRUE)
  rect(x_leg0, y_vals[1], x_leg1, y_vals[n_leg+1],
       col=NA, border=rgb(1,1,1,0.3), lwd=0.5, xpd=TRUE)
  at_days = c(0, round(max_day/3), round(2*max_day/3), max_day)
  at_y    = usr[3] + 0.1*(usr[4]-usr[3]) + (at_days/max_day)*0.8*(usr[4]-usr[3])
  for (k in seq_along(at_days))
    text(x_leg1 + 0.005*(usr[2]-usr[1]), at_y[k], labels=paste0("j.", at_days[k]),
         col=rgb(1,1,1,0.75), cex=0.45, adj=0, xpd=TRUE)
  text((x_leg0+x_leg1)/2, y_vals[n_leg+1] + 0.025*(usr[4]-usr[3]),
       labels="Jour\nde\ncollecte", col=rgb(1,1,1,0.6), cex=0.42, adj=0.5, xpd=TRUE, font=2)
  
  # Légende des symboles
  legend("bottomleft",
         legend   = c(sprintf("Cas récents (≤15j)  [n=%d]", sum(is_recent)),
                      sprintf("Cas anciens  [n=%d]", sum(!is_recent)),
                      "Zone infectée KDE 95%",
                      "Autoroute",
                      "Frontière"),
         col      = c(col_halo, cols_anim[30], col_kde, col_mway, col_border),
         pch      = c(16, 16, 15, NA, NA),
         lty      = c(NA, NA, NA,  1,  2),
         lwd      = c(NA, NA, NA,  2.0, 1.0),
         pt.cex   = c(1.4, 0.9, 1.5, NA, NA),
         cex      = 0.58,
         bty      = "n",
         text.col = rgb(1,1,1,0.85),
         bg       = rgb(0,0,0,0.45),
         box.col  = NA,
         inset    = c(0.01, 0.04))
  
  dev.off()
  cat(sprintf("  Frame %d/%d  (Jour %d)\r", i, length(days_seq), d))
}

cat("\nToutes les frames générées.\n")


# Compilation du GIF à partir des images générées

frame_files = sort(list.files(frames_dir, pattern="\\.png$", full.names=TRUE))
n_frames    = length(frame_files)

# Ajustement de la vitesse de l'animation pour un début et une fin plus lents
delays = rep(0.10, n_frames)
delays[1:min(10, n_frames)]            = 0.22
delays[max(1, n_frames-9):n_frames]    = 0.22
delays[max(1, n_frames-4):n_frames]    = 0.40

# Création du fichier final
gifski(png_files=frame_files, gif_file="animation_progression_cas.gif",
       width=1200, height=900, delay=delays, loop=TRUE)

cat("\nGIF sauvegardé : wavefront_progression_v2.gif\n")
cat(sprintf("  %d frames  |  Pas : 3j  |  1200x900 px\n", n_frames))

