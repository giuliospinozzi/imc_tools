# script.R

# Specifica il percorso completo del file args.txt
args_file <- "/home/rstudio/IMCDataAnalysis/args.txt"

# Leggi il file e ignora le righe che iniziano con '#'
args <- readLines(args_file)
args <- args[!grepl("^#", args)]  # Rimuovi le righe che iniziano con '#'
args <- args[args != ""]          # Rimuovi le righe vuote

# Verifica se hai ricevuto esattamente 8 argomenti
if (length(args) != 8) {
  stop("Il file deve contenere esattamente 8 linee non commentate: marker_list, ch_to_exclude, workingdir, scale_value, batch_correction, use_louvain_clusters, manual_clusters, e roi_list.")
}

# Assegna gli argomenti a variabili specifiche
marker_list <- strsplit(args[1], ",")[[1]]
ch_to_exclude <- args[2]
workingdir <- file.path("/home/rstudio/IMCDataAnalysis", args[3])
scale_value <- as.numeric(args[4])
batch_correction <- args[5]
use_louvain_clusters <- args[6] == "TRUE"
manual_clusters <- as.integer(args[7])
roi_list <- as.integer(strsplit(args[8], ",")[[1]])

# Verifica il valore di scale_value
if (is.na(scale_value) || scale_value < 1 || scale_value > 5) {
  stop("Il valore di scale_value deve essere un numero tra 1 e 5.")
}

# Stampa gli argomenti per verificare
print("Lista dei marker:")
print(marker_list)
print(paste("Canali da escludere:", ch_to_exclude))
print(paste("Directory di lavoro:", workingdir))
print(paste("Scale value:", scale_value))
print(paste("Batch correction:", batch_correction))
print(paste("Usare cluster di Louvain:", use_louvain_clusters))
print(paste("Numero di cluster manuale:", manual_clusters))
print(paste("Lista delle ROI:", roi_list))

# Crea la directory di lavoro
dir.create(workingdir, showWarnings = FALSE)
setwd(workingdir)
library(imcRtools)
library(cytomapper)
library(viridis)
library(dittoSeq)
library(tidyverse)
library(ggrepel)
library(EBImage)
library(CATALYST)
library(BiocParallel)

####################MODIFICHE PER COMPENSAZIONE##################################
path <- "/home/rstudio/IMCDataAnalysis"
spe <- read_steinbock(path)
colnames(spe) <- paste0(spe$sample_id, "_", spe$ObjectNumber) #Saves ROI id sa sample_id

# Estrai i nomi delle colonne e crea una singola stringa separata da '|'
colnames_spe <- colnames(spe)
unique_colnames <- unique(sub("_\\d+$", "", colnames_spe))
regexp <- paste(unique_colnames, collapse = "|")

# Associa il patient_id
spe$patient_id <- str_extract(spe$sample_id, regexp)

# Stampa i risultati per verificare
print(table(spe$patient_id))

# Salva i risultati in un file
write.table(table(spe$patient_id), "n_cells.tsv", sep = ",", col.names = NA)

#indication <- meta$Indication[match(spe$patient_id, meta$`Sample ID`)]
rowData(spe)$use_channel <- !grepl(ch_to_exclude, rownames(spe))
## Performs scaling using arcsin function and a given scaling value and re-displays data
scale_value = 2
assay(spe, "exprs") <- asinh(counts(spe)/scale_value)
masks <- loadImages(paste(path,'/masks/',sep=''), as.is = TRUE) #edited as.is as Dario suggested 
images <- loadImages(paste(path,'/img/',sep=''))
channelNames(images) <- rownames(spe)

##########################Questa cosa non c'è sul nature protocol#############

all.equal(names(images), names(masks))

##############################################################################
patient_id <- str_extract(names(images), regexp)

#indication <- meta$Indication[match(patient_id, meta$`Sample ID`)]

#EDITATO COME SUL NATURE PROTOCOL
mcols(images) <- mcols(masks) <- DataFrame(sample_id = names(images),
                                           patient_id = patient_id)

#spillover correction

sce <- readSCEfromTXT("/home/rstudio/IMCDataAnalysis/compensation/") 
assay(sce, "exprs") <- asinh(counts(sce)/5)
plotSpotHeatmap(sce)
sce2 <- binAcrossPixels(sce, bin_size = 10)
library(CATALYST)
bc_key <- as.numeric(unique(sce$sample_mass)) 
bc_key <- bc_key[order(bc_key)]
sce <- assignPrelim(sce, bc_key = bc_key)
sce <- estCutoffs(sce)
sce <- applyCutoffs(sce)
library(pheatmap)
cur_table <- table(sce$bc_id, sce$sample_mass)
# Visualize the correctly and incorrectly assigned pixels 
pheatmap(log10(cur_table + 1),
         cluster_rows = FALSE,
         cluster_cols = FALSE)
# Compute the fraction of unassigned pixels per spot
cur_table["0",] / colSums(cur_table)
# Filter pixels that were assigned to a mass other than the spotted mass 
sce <- filterPixels(sce, minevents = 40, correct_pixels = TRUE)
sce <- computeSpillmat(sce)
sm <- metadata(sce)$spillover_matrix
library(dittoSeq) 
library(patchwork)
# Specify the channel_name entry for use with CATALYST 
rowData(spe)$channel_name <- paste0(rowData(spe)$channel, "Di")
# Add the missing metal isotope to the 
isotope_list <- CATALYST::isotope_list 
isotope_list$Ar <- 80

spe <- compCytof(spe, sm, transform = TRUE, cofactor = 1, isotope_list = isotope_list, overwrite = FALSE)

# Replace uncompensated assays
assay(spe, "counts") <- assay(spe, "compcounts")
assay(spe, "exprs") <- assay(spe, "compexprs")
assay(spe, "compcounts") <- assay(spe, "compexprs") <- NULL
# Use mass tags as channel names
channelNames(images) <- rowData(spe)$channel_name
# Adapt spillover matrix to retain channels contained # in the multi-channel images
adapted_sm <- adaptSpillmat(sm, channelNames(images),isotope_list = isotope_list)
# Perform image compensation
channelNames(images) <- gsub("ArAr80Di", "Ar80Di", channelNames(images))
colnames(adapted_sm) <- gsub("ArAr80Di", "Ar80Di", colnames(adapted_sm))
################################# ERRORE #######################################

colonne_da_tenere <- intersect(channelNames(images), colnames(adapted_sm))
adapted_sm <- adapted_sm[, colonne_da_tenere]

images_comp <- compImage(images, adapted_sm, BPPARAM =  MulticoreParam())

##############################ESTRAZIONE IMMAGINI AUTOMATIZZATA##################################################
#SELEZIONA PRIMA IL NUMERO DELLA ROI DA SALVARE, MODIFICA I NOMI AGGIUNGENDO IL NUMERO DELLA ROI, POI ESEGUI

# Ottieni i nomi dei canali
channels <- channelNames(images)

# Itera su ciascun canale e ROI per generare i plot per before e after
for (roi in roi_list) {
  for (channel in channels) {
    # Tentativo di generare il plot prima della compensazione
    tryCatch({
      png(filename = paste0("plot_ROI", roi, "_", channel, "_before.png"), width = 800, height = 600)
      plotPixels(
        images[roi], 
        colour_by = channel, 
        image_title = list(text = paste(channel, "- before"), position = "topleft"),
        legend = NULL, 
        bcg = setNames(list(c(0, 4, 2)), channel)  # Imposta il nome corretto del canale
      )
      dev.off()
    }, error = function(e) {
      message(paste("Errore nel plotting prima della compensazione per ROI", roi, "e canale", channel, ":", e$message))
      dev.off()
    })
    
    # Tentativo di generare il plot dopo la compensazione
    tryCatch({
      png(filename = paste0("plot_ROI", roi, "_", channel, "_after.png"), width = 800, height = 600)
      plotPixels(
        images_comp[roi], 
        colour_by = channel, 
        image_title = list(text = paste(channel, "- after"), position = "topleft"),
        legend = NULL, 
        bcg = setNames(list(c(0, 4, 2)), channel)  # Imposta il nome corretto del canale
      )
      dev.off()
    }, error = function(e) {
      message(paste("Errore nel plotting dopo la compensazione per ROI", roi, "e canale", channel, ":", e$message))
      dev.off()
    })
  }
}

# Switch back to using target names as channel names 
channelNames(images_comp) <- rownames(spe)
set.seed(20220118)
img_ids <- sample(seq_len(length(images_comp)), 3)
cur_images <- images_comp[img_ids]
# Normalize each channel between 0 and 1
cur_images <- cytomapper::normalize(cur_images, separateImages = TRUE) 
# Clip channel intensities at 0 and 0.2



dittoPlot(spe, var = "area", group.by = "patient_id", plots = "boxplot") + ylab("Cell area") + xlab("")
spe <- spe[,spe$area >= 5]
dittoPlot(spe, var = "area", group.by = "sample_id", plots = "boxplot") + ylab("Cell area") + xlab("")
spe <- spe[,spe$area >= 5]







# Itera attraverso la lista di marker
for (marker_to_use in marker_list) {
  # Calcola la trasformazione
  assay(spe, "exprs") <- asinh(counts(spe) / scale_value)
  
  # Genera il plot
  plot <- dittoRidgePlot(spe, var = marker_to_use, group.by = "patient_id", assay = "exprs") +
    ggtitle(paste(marker_to_use, " - post transformation", sep = ''))
  
  # Salva il plot
  ggsave(filename = paste0(marker_to_use, "_plot.pdf"), plot = plot,
         width = 10,  # Imposta la larghezza desiderata in pollici
         height = 8)  # Imposta l'altezza desiderata in pollici)
  
}


## Saves data
#saveRDS(spe, paste(path,"/spe.rds",sep=''))
#saveRDS(images, paste(path,"/images.rds",sep=''))
#saveRDS(masks, paste(path,"/masks.rds",sep=''))

## Heatmap of  all the good working markers for each cell - overlayed with ROI they belong to (sampled to just 2000 cells for speed)
cells_to_sample = 2000

cur_cells <- sample(seq_len(ncol(spe)), cells_to_sample)

library(colorspace)
A<-c((0.18995),(0.19483),(0.19956),(0.20415),(0.2086),(0.21291),(0.21708),(0.22111),(0.225),(0.22875),(0.23236),(0.23582),(0.23915),(0.24234),(0.24539),(0.2483),(0.25107),(0.25369),(0.25618),(0.25853),(0.26074),(0.2628),(0.26473),(0.26652),(0.26816),(0.26967),(0.27103),(0.27226),(0.27334),(0.27429),(0.27509),(0.27576),(0.27628),(0.27667),(0.27691),(0.27701),(0.27698),(0.2768),(0.27648),(0.27603),(0.27543),(0.27469),(0.27381),(0.27273),(0.27106),(0.26878),(0.26592),(0.26252),(0.25862),(0.25425),(0.24946),(0.24427),(0.23874),(0.23288),(0.22676),(0.22039),(0.21382),(0.20708),(0.20021),(0.19326),(0.18625),(0.17923),(0.17223),(0.16529),(0.15844),(0.15173),(0.14519),(0.13886),(0.13278),(0.12698),(0.12151),(0.11639),(0.11167),(0.10738),(0.10357),(0.10026),(0.0975),(0.09532),(0.09377),(0.09287),(0.09267),(0.0932),(0.09451),(0.09662),(0.09958),(0.10342),(0.10815),(0.11374),(0.12014),(0.12733),(0.13526),(0.14391),(0.15323),(0.16319),(0.17377),(0.18491),(0.19659),(0.20877),(0.22142),(0.23449),(0.24797),(0.2618),(0.27597),(0.29042),(0.30513),(0.32006),(0.33517),(0.35043),(0.36581),(0.38127),(0.39678),(0.41229),(0.42778),(0.44321),(0.45854),(0.47375),(0.48879),(0.50362),(0.51822),(0.53255),(0.54658),(0.56026),(0.57357),(0.58646),(0.59891),(0.61088),(0.62233),(0.63323),(0.64362),(0.65394),(0.66428),(0.67462),(0.68494),(0.69525),(0.70553),(0.71577),(0.72596),(0.7361),(0.74617),(0.75617),(0.76608),(0.77591),(0.78563),(0.79524),(0.80473),(0.8141),(0.82333),(0.83241),(0.84133),(0.8501),(0.85868),(0.86709),(0.8753),(0.88331),(0.89112),(0.8987),(0.90605),(0.91317),(0.92004),(0.92666),(0.93301),(0.93909),(0.94489),(0.95039),(0.9556),(0.96049),(0.96507),(0.96931),(0.97323),(0.97679),(0.98),(0.98289),(0.98549),(0.98781),(0.98986),(0.99163),(0.99314),(0.99438),(0.99535),(0.99607),(0.99654),(0.99675),(0.99672),(0.99644),(0.99593),(0.99517),(0.99419),(0.99297),(0.99153),(0.98987),(0.98799),(0.9859),(0.9836),(0.98108),(0.97837),(0.97545),(0.97234),(0.96904),(0.96555),(0.96187),(0.95801),(0.95398),(0.94977),(0.94538),(0.94084),(0.93612),(0.93125),(0.92623),(0.92105),(0.91572),(0.91024),(0.90463),(0.89888),(0.89298),(0.88691),(0.88066),(0.87422),(0.8676),(0.86079),(0.8538),(0.84662),(0.83926),(0.83172),(0.82399),(0.81608),(0.80799),(0.79971),(0.79125),(0.7826),(0.77377),(0.76476),(0.75556),(0.74617),(0.73661),(0.72686),(0.71692),(0.7068),(0.6965),(0.68602),(0.67535),(0.66449),(0.65345),(0.64223),(0.63082),(0.61923),(0.60746),(0.5955),(0.58336),(0.57103),(0.55852),(0.54583),(0.53295),(0.51989),(0.50664),(0.49321),(0.4796))
B<-c((0.07176),(0.08339),(0.09498),(0.10652),(0.11802),(0.12947),(0.14087),(0.15223),(0.16354),(0.17481),(0.18603),(0.1972),(0.20833),(0.21941),(0.23044),(0.24143),(0.25237),(0.26327),(0.27412),(0.28492),(0.29568),(0.30639),(0.31706),(0.32768),(0.33825),(0.34878),(0.35926),(0.3697),(0.38008),(0.39043),(0.40072),(0.41097),(0.42118),(0.43134),(0.44145),(0.45152),(0.46153),(0.47151),(0.48144),(0.49132),(0.50115),(0.51094),(0.52069),(0.5304),(0.54015),(0.54995),(0.55979),(0.56967),(0.57958),(0.5895),(0.59943),(0.60937),(0.61931),(0.62923),(0.63913),(0.64901),(0.65886),(0.66866),(0.67842),(0.68812),(0.69775),(0.70732),(0.7168),(0.7262),(0.73551),(0.74472),(0.75381),(0.76279),(0.77165),(0.78037),(0.78896),(0.7974),(0.80569),(0.81381),(0.82177),(0.82955),(0.83714),(0.84455),(0.85175),(0.85875),(0.86554),(0.87211),(0.87844),(0.88454),(0.8904),(0.896),(0.90142),(0.90673),(0.91193),(0.91701),(0.92197),(0.9268),(0.93151),(0.93609),(0.94053),(0.94484),(0.94901),(0.95304),(0.95692),(0.96065),(0.96423),(0.96765),(0.97092),(0.97403),(0.97697),(0.97974),(0.98234),(0.98477),(0.98702),(0.98909),(0.99098),(0.99268),(0.99419),(0.99551),(0.99663),(0.99755),(0.99828),(0.99879),(0.9991),(0.99919),(0.99907),(0.99873),(0.99817),(0.99739),(0.99638),(0.99514),(0.99366),(0.99195),(0.98999),(0.98775),(0.98524),(0.98246),(0.97941),(0.9761),(0.97255),(0.96875),(0.9647),(0.96043),(0.95593),(0.95121),(0.94627),(0.94113),(0.93579),(0.93025),(0.92452),(0.91861),(0.91253),(0.90627),(0.89986),(0.89328),(0.88655),(0.87968),(0.87267),(0.86553),(0.85826),(0.85087),(0.84337),(0.83576),(0.82806),(0.82025),(0.81236),(0.80439),(0.79634),(0.78823),(0.78005),(0.77181),(0.76352),(0.75519),(0.74682),(0.73842),(0.73),(0.7214),(0.7125),(0.7033),(0.69382),(0.68408),(0.67408),(0.66386),(0.65341),(0.64277),(0.63193),(0.62093),(0.60977),(0.59846),(0.58703),(0.57549),(0.56386),(0.55214),(0.54036),(0.52854),(0.51667),(0.50479),(0.49291),(0.48104),(0.4692),(0.4574),(0.44565),(0.43399),(0.42241),(0.41093),(0.39958),(0.38836),(0.37729),(0.36638),(0.35566),(0.34513),(0.33482),(0.32473),(0.31489),(0.3053),(0.29599),(0.28696),(0.27824),(0.26981),(0.26152),(0.25334),(0.24526),(0.2373),(0.22945),(0.2217),(0.21407),(0.20654),(0.19912),(0.19182),(0.18462),(0.17753),(0.17055),(0.16368),(0.15693),(0.15028),(0.14374),(0.13731),(0.13098),(0.12477),(0.11867),(0.11268),(0.1068),(0.10102),(0.09536),(0.0898),(0.08436),(0.07902),(0.0738),(0.06868),(0.06367),(0.05878),(0.05399),(0.04931),(0.04474),(0.04028),(0.03593),(0.03169),(0.02756),(0.02354),(0.01963),(0.01583))
C<-c((0.23217),(0.26149),(0.29024),(0.31844),(0.34607),(0.37314),(0.39964),(0.42558),(0.45096),(0.47578),(0.50004),(0.52373),(0.54686),(0.56942),(0.59142),(0.61286),(0.63374),(0.65406),(0.67381),(0.693),(0.71162),(0.72968),(0.74718),(0.76412),(0.7805),(0.79631),(0.81156),(0.82624),(0.84037),(0.85393),(0.86692),(0.87936),(0.89123),(0.90254),(0.91328),(0.92347),(0.93309),(0.94214),(0.95064),(0.95857),(0.96594),(0.97275),(0.97899),(0.98461),(0.9893),(0.99303),(0.99583),(0.99773),(0.99876),(0.99896),(0.99835),(0.99697),(0.99485),(0.99202),(0.98851),(0.98436),(0.97959),(0.97423),(0.96833),(0.9619),(0.95498),(0.94761),(0.93981),(0.93161),(0.92305),(0.91416),(0.90496),(0.8955),(0.8858),(0.8759),(0.86581),(0.85559),(0.84525),(0.83484),(0.82437),(0.81389),(0.80342),(0.79299),(0.78264),(0.7724),(0.7623),(0.75237),(0.74265),(0.73316),(0.72393),(0.715),(0.70599),(0.69651),(0.6866),(0.67627),(0.66556),(0.65448),(0.64308),(0.63137),(0.61938),(0.60713),(0.59466),(0.58199),(0.56914),(0.55614),(0.54303),(0.52981),(0.51653),(0.50321),(0.48987),(0.47654),(0.46325),(0.45002),(0.43688),(0.42386),(0.41098),(0.39826),(0.38575),(0.37345),(0.3614),(0.34963),(0.33816),(0.32701),(0.31622),(0.30581),(0.29581),(0.28623),(0.27712),(0.26849),(0.26038),(0.2528),(0.24579),(0.23937),(0.23356),(0.22835),(0.2237),(0.2196),(0.21602),(0.21294),(0.21032),(0.20815),(0.2064),(0.20504),(0.20406),(0.20343),(0.20311),(0.2031),(0.20336),(0.20386),(0.20459),(0.20552),(0.20663),(0.20788),(0.20926),(0.21074),(0.2123),(0.21391),(0.21555),(0.21719),(0.2188),(0.22038),(0.22188),(0.22328),(0.22456),(0.2257),(0.22667),(0.22744),(0.228),(0.22831),(0.22836),(0.22811),(0.22754),(0.22663),(0.22536),(0.22369),(0.22161),(0.21918),(0.2165),(0.21358),(0.21043),(0.20706),(0.20348),(0.19971),(0.19577),(0.19165),(0.18738),(0.18297),(0.17842),(0.17376),(0.16899),(0.16412),(0.15918),(0.15417),(0.1491),(0.14398),(0.13883),(0.13367),(0.12849),(0.12332),(0.11817),(0.11305),(0.10797),(0.10294),(0.09798),(0.0931),(0.08831),(0.08362),(0.07905),(0.07461),(0.07031),(0.06616),(0.06218),(0.05837),(0.05475),(0.05134),(0.04814),(0.04516),(0.04243),(0.03993),(0.03753),(0.03521),(0.03297),(0.03082),(0.02875),(0.02677),(0.02487),(0.02305),(0.02131),(0.01966),(0.01809),(0.0166),(0.0152),(0.01387),(0.01264),(0.01148),(0.01041),(0.00942),(0.00851),(0.00769),(0.00695),(0.00629),(0.00571),(0.00522),(0.00481),(0.00449),(0.00424),(0.00408),(0.00401),(0.00401),(0.0041),(0.00427),(0.00453),(0.00486),(0.00529),(0.00579),(0.00638),(0.00705),(0.0078),(0.00863),(0.00955),(0.01055))
turbo_colormap_data<-cbind(A,B,C) 
turbo_colormap_data_sRGB<-sRGB(turbo_colormap_data)
turbo_colormap_data_HEX = hex(turbo_colormap_data_sRGB)

col_sample_id <- setNames(dittoColors(2)[1:length(unique(spe$patient_id))],
                          unique(spe$patient_id))


# p<-dittoHeatmap(spe[,cur_cells], treeheight_row=0,  treeheight_col=0, genes = rownames(spe)[rowData(spe)$use_channel],
#              assay = "exprs", cluster_cols = F, scale = "none",
#              heatmap.colors = colorRampPalette(colors = turbo_colormap_data_HEX)(100), annot.by = "sample_id",
#              annotation_colors = list(sample_id = col_sample_id))

callback_rev <- function(hc, mat){
  hc$order <- sort(hc$order)
  hc
}

p<-dittoHeatmap(spe[,cur_cells], treeheight_row=0,  treeheight_col=0, genes = rownames(spe)[rowData(spe)$use_channel],
                assay = "exprs", cluster_cols = F, scale = "none",
                heatmap.colors = colorRampPalette(colors = turbo_colormap_data_HEX)(100), annot.by = "patient_id",
                annotation_colors = list(patient_id = col_sample_id),clustering_callback = callback_rev)


#p<-dittoHeatmap(spe[,cur_cells], genes = rownames(spe)[rowData(spe)$use_channel],
#             assay = "exprs", cluster_cols = TRUE, scale = "none",
#             heatmap.colors = colorRampPalette(colors = turbo_colormap_data_HEX)(100), annot.by = "sample_id")
ggsave(filename="Heatmap-All-Markers.pdf", p, width=10, height=12)
ggsave(filename="Heatmap-All-Markers.png", p, width=10, height=12)

p<-dittoHeatmap(spe[,cur_cells], treeheight_row=50,  treeheight_col=50, genes = rownames(spe)[rowData(spe)$use_channel],
                assay = "exprs", cluster_cols = F, scale = "none",
                heatmap.colors = colorRampPalette(colors = turbo_colormap_data_HEX)(100), annot.by = "patient_id",
                annotation_colors = list(patient_id = col_sample_id))


#p<-dittoHeatmap(spe[,cur_cells], genes = rownames(spe)[rowData(spe)$use_channel],
#             assay = "exprs", cluster_cols = TRUE, scale = "none",
#             heatmap.colors = colorRampPalette(colors = turbo_colormap_data_HEX)(100), annot.by = "sample_id")
ggsave(filename="Heatmap-All-Markers-Dendrogramma.pdf", p, width=10, height=12)
ggsave(filename="Heatmap-All-Markers-Dendrogramma.png", p, width=10, height=12)

## Plots Signal intensity vs. SNR for all channels (useful to see if any did not work)
library(ggrepel)
cur_snr <- lapply(images, function(img){
  mat <- apply(img, 3, function(ch){
    # Otsu threshold
    thres <- otsu(ch, range = c(min(ch), max(ch)))
    # Signal-to-noise ratio
    snr <- mean(ch[ch > thres]) / mean(ch[ch <= thres])
    # Signal intensity
    ps <- mean(ch[ch > thres])
    
    return(c(snr = snr, ps = ps))
  })
  t(mat) %>% as.data.frame() %>% 
    mutate(marker = colnames(mat)) %>% 
    pivot_longer(cols = c(snr, ps))
})

cur_snr <- do.call(rbind, cur_snr)
cur_snr %>% 
  group_by(marker, name) %>%
  summarize(mean = mean(value),
            ci = qnorm(0.975)*sd(value)/sqrt(n())) %>%
  pivot_wider(names_from = name, values_from = c(mean, ci)) %>%
  ggplot() +
  geom_point(aes(log2(mean_ps), log2(mean_snr))) +
  geom_label_repel(aes(log2(mean_ps), log2(mean_snr), label = marker)) +
  theme_minimal(base_size = 15) + ylab("Signal-to-noise ratio [log2]") +
  xlab("Signal intensity [log2]")
p <- cur_snr %>% 
  group_by(marker, name) %>%
  summarize(mean = mean(value),
            ci = qnorm(0.975)*sd(value)/sqrt(n())) %>%
  pivot_wider(names_from = name, values_from = c(mean, ci)) %>%
  ggplot() +
  geom_point(aes(log2(mean_ps), log2(mean_snr))) +
  geom_label_repel(aes(log2(mean_ps), log2(mean_snr), label = marker)) +
  theme_minimal(base_size = 15) + ylab("Signal-to-noise ratio [log2]") +
  xlab("Signal intensity [log2]")

# Salva il plot
ggsave(filename = paste0("Signal_intensity.pdf"), plot = p,width=10, height=10)
## For each ROI, plots how much of the area is covered by cells (good indication of cell density)
colData(spe) %>%
  as.data.frame() %>%
  group_by(patient_id) %>%
  summarize(cell_area = sum(area),
            no_pixels = mean(width_px) * mean(height_px)) %>%
  mutate(covered_area = cell_area / no_pixels) %>%
  ggplot() +
  geom_point(aes(reorder(patient_id,covered_area), covered_area)) + 
  theme_minimal(base_size = 15) +
  ylim(c(0, 1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("% covered area") + xlab("")

## Plots average cell area per ROI (useful to see if there are large-scale ROI differences)
colData(spe) %>%
  as.data.frame() %>%
  group_by(patient_id) %>%
  ggplot() +
  geom_boxplot(aes(patient_id, area)) +
  theme_minimal(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("Cell area") + xlab("")
p <- colData(spe) %>%
  as.data.frame() %>%
  group_by(patient_id) %>%
  ggplot() +
  geom_boxplot(aes(patient_id, area)) +
  theme_minimal(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("Cell area") + xlab("")
# Salva il plot
ggsave(filename = paste0("Cell_area.pdf"), plot = p)
summary(spe$area)


## Plots average cell density per ROI (useful to seeimages## Plots average cell density per ROI (useful to see if there are large-scale ROI differences)
colData(spe) %>%
  as.data.frame() %>%
  group_by(patient_id) %>%
  summarize(cell_count = n(),
            no_pixels = mean(width_px) * mean(height_px)) %>%
  mutate(cells_per_mm2 = cell_count/(no_pixels/1000000)) %>%
  ggplot() +
  geom_point(aes(patient_id, cells_per_mm2)) + 
  theme_minimal(base_size = 15)  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("Cells per mm2") + xlab("")
p<- colData(spe) %>%
  as.data.frame() %>%
  group_by(patient_id) %>%
  summarize(cell_count = n(),
            no_pixels = mean(width_px) * mean(height_px)) %>%
  mutate(cells_per_mm2 = cell_count/(no_pixels/1000000)) %>%
  ggplot() +
  geom_point(aes(patient_id, cells_per_mm2)) + 
  theme_minimal(base_size = 15)  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("Cells per mm2") + xlab("")
ggsave(filename = paste0("Cell_per_mm2.pdf"), plot = p)
sum(spe$area < 25)
spe <- spe[,spe$area >= 25]


## Plots ridge plots for all used channels (using rescaled expression)
plot <- multi_dittoPlot(spe, vars = rownames(spe)[rowData(spe)$use_channel],
                        group.by = "patient_id", plots = c("ridgeplot"), 
                        assay = "exprs")
# Salva il plot
ggsave(filename = paste0("Distribuzioni", "_multiplot.pdf"), plot = plot,
       width = 10,  # Imposta la larghezza desiderata in pollici
       height = 20)  # Imposta l'altezza desiderata in pollici


## Runs UMAP dimensionality reduction using all the channels defined as working, and the re-scaled expression
library(scater)
set.seed(220225)
spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs") 

## Plots ROI on UMAP
library(patchwork)

# visualize sample id 
p<-dittoDimPlot(spe, var = "patient_id", reduction.use = "UMAP", size = 0.4) + 
  ggtitle("ROI ID on UMAP")
dittoDimPlot(spe, var = "patient_id", reduction.use = "UMAP", size = 0.4) + 
  ggtitle("ROI ID on UMAP")
ggsave(filename="ROI-ID-UMAP.pdf", p, width=10, height=6)
ggsave(filename="ROI-ID-UMAP.png", p, width=10, height=6)


## Visualize marker expression for any given marker on UMAP
#use_marker = "Caspase3"

for (m in marker_list) {
  p<-dittoDimPlot(spe, var = m, reduction.use = "UMAP", 
                  assay = "exprs", size = 0.2) + 
    scale_color_viridis(option = "H",name=m) +
    ggtitle(paste(m," expression on UMAP",sep=''))
  ggsave(filename=paste0("UMAP_",m,".pdf"), p)
}

####################### PROVA AUTOMAZIONE CLUSTER #######################


library(batchelor)
library(scater)
library(bluster)
library(BiocParallel)
library(dittoSeq)
library(viridis)

set.seed(220228)

#########################################################################

# Controlla se deve essere applicata la correzione batch
if (batch_correction == "Si") {
  # Correzione batch con fastMNN
  out <- fastMNN(spe, batch = spe$patient_id,
                 auto.merge = TRUE,
                 subset.row = rowData(spe)$use_channel,
                 assay.type = "exprs")
  
  # Memorizza le embedding corrette nell'oggetto SPE
  reducedDim(spe, "fastMNN") <- reducedDim(out, "corrected")
  
  # Esegui UMAP sulle embedding corrette
  spe <- runUMAP(spe, dimred = "fastMNN", name = "UMAP_mnnCorrected")
  
  # Visualizza UMAP corretto
  p_umap_corrected <- dittoDimPlot(spe, var = "patient_id",
                                   reduction.use = "UMAP_mnnCorrected", size = 0.2) +
    ggtitle("Patient ID on UMAP after correction")
  ggsave(filename = "UMAP_mnnCorrected.pdf", p_umap_corrected, width = 10, height = 6)
  ggsave(filename = "UMAP_mnnCorrected.png", p_umap_corrected, width = 10, height = 6)
  
  # Seleziona le embedding corrette per il clustering
  mat <- reducedDim(spe, "fastMNN")
} else {
  # Esegui UMAP sui dati originali senza correzione batch
  spe <- runUMAP(spe, exprs_values = "exprs", name = "UMAP_original")
  
  # Visualizza UMAP senza correzione
  p_umap_original <- dittoDimPlot(spe, var = "patient_id",
                                  reduction.use = "UMAP_original", size = 0.2) +
    ggtitle("Patient ID on UMAP without correction")
  ggsave(filename = "UMAP_original.pdf", p_umap_original, width = 10, height = 6)
  ggsave(filename = "UMAP_original.png", p_umap_original, width = 10, height = 6)
  
  # Seleziona le embedding UMAP per il clustering
  mat <- reducedDim(spe, "UMAP_original")
}

# Esegui il cluster sweep con Louvain
combinations <- clusterSweep(mat, BLUSPARAM = SNNGraphParam(),
                             k = c(10L, 20L),
                             type = c("rank", "jaccard"),
                             cluster.fun = "louvain",
                             BPPARAM = SerialParam(RNGseed = 230214))

# Calcola la larghezza media della silhouette per ogni combinazione di parametri
sil <- vapply(as.list(combinations$clusters),
              function(x) mean(approxSilhouette(mat, x)$width), 0)

# Trova la combinazione di parametri ottimale
best_index <- which.max(sil)
best_combination <- combinations$clusters[[best_index]]

# Assegna i cluster ottenuti da Louvain
spe$nn_clusters <- best_combination

# Conta il numero di cluster ottenuti
num_clusters <- length(unique(best_combination))

# Genera etichette numeriche crescenti per ogni cluster
cluster_labels <- as.character(seq_len(num_clusters))
names(cluster_labels) <- as.character(seq_len(num_clusters))

# Assegna le etichette ai cluster
cluster_celltype <- recode(spe$nn_clusters, !!!setNames(cluster_labels, cluster_labels))
spe$cluster_celltype <- cluster_celltype

# Visualizza la heatmap dei cluster Louvain
cur_cells <- sample(seq_len(ncol(spe)), 2000)
p <- dittoHeatmap(spe[,cur_cells],
                  genes = rownames(spe)[rowData(spe)$use_channel],
                  assay = "exprs", scale = "none",
                  heatmap.colors = viridis(100),
                  annot.by = c("nn_clusters", "patient_id"))

ggsave(filename = "Heatmap-Louvain-Clusters.pdf", p, width = 10, height = 14)
ggsave(filename = "Heatmap-Louvain-Clusters.png", p, width = 10, height = 14)

# Visualizza i cluster Louvain sull'UMAP corretto o non corretto
umap_name <- if(batch_correction == "Si") "UMAP_mnnCorrected" else "UMAP_original"

p2 <- dittoDimPlot(spe, var = "cluster_celltype", reduction.use = umap_name, size = 0.4, do.label = TRUE) +
  ggtitle("Cluster cell types on UMAP")

# Visualizza il plot
print(p2)

# Salva il plot
ggsave(filename = "Cluster_cell_types_on_UMAP.pdf", p2, width = 10, height = 6)

#########################################################################

# Esegui il clustering con FlowSOM e ConsensusClusterPlus, utilizzando i cluster Louvain
spe <- cluster(spe, 
               features = rownames(spe)[rowData(spe)$use_channel],
               maxK = 45,  # Numero massimo di cluster da considerare
               seed = 220410)  # Seed per garantire la riproducibilità

# Verifica se i risultati del clustering sono stati generati correttamente
if (is.null(metadata(spe)$cluster_codes)) {
  stop("Errore: Il clustering non è stato eseguito correttamente. Assicurati che la funzione cluster abbia funzionato.")
}

# Usa il numero di cluster trovato con Louvain
clusters_to_show <- num_clusters

# Verifica se esiste un cluster con il nome specificato
cluster_name <- paste('meta', clusters_to_show, sep = '')
if (!(cluster_name %in% names(metadata(spe)$cluster_codes))) {
  stop(paste("Errore: Il cluster specificato non esiste:", cluster_name))
}

# Assegna i cluster trovati come cluster SOM
spe$som_clusters <- cluster_ids(spe, cluster_name)

# Visualizza i cluster SOM sull'UMAP
p <- dittoDimPlot(spe, var = "som_clusters", 
                  reduction.use = "UMAP", size = 0.4,
                  do.label = TRUE) +
  ggtitle("SOM clusters expression on UMAP")
ggsave(filename = "Clusters-UMAP.pdf", p, width = 10, height = 6)
ggsave(filename = "Clusters-UMAP.png", p, width = 10, height = 6)

# Visualizza la heatmap dei cluster SOM
p <- dittoHeatmap(spe, 
                  genes = rownames(spe)[rowData(spe)$use_channel],
                  assay = "exprs", scale = "none", treeheight_row = 0,
                  heatmap.colors = colorRampPalette(colors = turbo_colormap_data_HEX)(100), 
                  annot.by = c("som_clusters", "patient_id"), 
                  annotation_colors = list(patient_id = col_sample_id), 
                  clustering_callback = callback_rev)

ggsave(filename = "Heatmap-Som-Clusters.pdf", p, width = 10, height = 14)
ggsave(filename = "Heatmap-Som-Clusters.png", p, width = 10, height = 14)
remove(p)

# Visualizza la heatmap dei cluster SOM con dendrogramma
p <- dittoHeatmap(spe, 
                  genes = rownames(spe)[rowData(spe)$use_channel],
                  assay = "exprs", scale = "none", treeheight_row = 50, treeheight_col = 50,
                  heatmap.colors = colorRampPalette(colors = turbo_colormap_data_HEX)(100), 
                  annot.by = c("som_clusters", "patient_id"))

ggsave(filename = "Heatmap-Som-Clusters-wDendrogram.pdf", p, width = 10, height = 14)
ggsave(filename = "Heatmap-Som-Clusters-wDendrogram.png", p, width = 10, height = 14)
remove(p)

#Vediamo
library(scuttle)

## aggregate by cell type
celltype_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),  
                                      ids = spe$som_clusters, 
                                      statistics = "mean",
                                      use.assay.type = "exprs", 
                                      subset.row = rownames(spe)[rowData(spe)$marker_class == "type"])

# No scaling
p<- dittoHeatmap(celltype_mean,
                 assay = "exprs", 
                 cluster_cols = TRUE, 
                 scale = "none",
                 heatmap.colors = viridis(100),
                 annot.by = c("som_clusters", "ncells"),
                 annotation_colors = list(celltype = metadata(spe)$color_vectors$celltype,
                                          ncells = plasma(100)))

ggsave(filename="Heatmap_no_scaling-UMAP.pdf", p, width=10, height=10)


# Visualize the results (heatmap)
library(viridis)
set.seed(220619)
cur_cells <- sample(seq_len(ncol(spe)), 2000) 
p <- dittoHeatmap(spe[, cur_cells], 
                  genes = rownames(spe)[rowData(spe)$use_channel], 
                  assay = "exprs", 
                  scale = "none",
                  heatmap.colors = viridis(100),
                  annot.by = c("nn_clusters", "patient_id"))
ggsave(filename="Heatmap_Clusters-UMAP.pdf", p, width=10, height=14)
#vediamo
library(scuttle)

## aggregate by cell type
celltype_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),  
                                      ids = spe$nn_clusters, 
                                      statistics = "mean",
                                      use.assay.type = "exprs", 
                                      subset.row = rownames(spe)[rowData(spe)$use_channel])

# No scaling
p <- dittoHeatmap(celltype_mean,
                  assay = "exprs", 
                  cluster_cols = TRUE, 
                  scale = "none",
                  heatmap.colors = viridis(100),
                  annot.by = c("nn_clusters", "ncells"),
                  annotation_colors = list(celltype = metadata(spe)$color_vectors$celltype,
                                           ncells = plasma(100)))
ggsave(filename="Heatmap_Clusters-no-scaling-UMAP.pdf", p, width=10, height=14)
# Scaled to max
p<- dittoHeatmap(celltype_mean,
                 assay = "exprs", 
                 cluster_cols = TRUE, 
                 scaled.to.max = TRUE,
                 heatmap.colors.max.scaled = viridis(100),
                 annot.by = c("nn_clusters", "ncells"),
                 annotation_colors = list(celltype = metadata(spe)$color_vectors$celltype,
                                          ncells = plasma(100)))
ggsave(filename="Heatmap_Clusters-scaled-to-max-UMAP.pdf", p, width=10, height=14)
# Z score scaled
p<-dittoHeatmap(celltype_mean,
                assay = "exprs", 
                cluster_cols = TRUE, 
                annot.by = c("nn_clusters", "ncells"),
                annotation_colors = list(celltype = metadata(spe)$color_vectors$celltype,
                                         ncells = viridis(100)))
ggsave(filename="Heatmap_Clusters_z-score-UMAP.pdf", p, width=10, height=14)




library(ggplot2)


# Numero di cluster da usare per FlowSOM e Steinbock graph
if (use_louvain_clusters) {
  clusters_to_use <- spe$nn_clusters  # Usa i cluster trovati da Louvain
} else {
  clusters_to_use <- cluster_ids(spe, paste('meta', manual_clusters, sep = ''))  # Usa il numero di cluster specificato manualmente
}

# Assegna i cluster SOM direttamente dai cluster di Louvain (o manuali)
spe$som_clusters <- clusters_to_use

# Genera una palette di colori sufficiente per coprire tutti i cluster SOM
num_clusters_total <- length(unique(spe$som_clusters))
color_palette <- scales::hue_pal()(num_clusters_total)

# Associa i colori generati ai cluster
cluster_colors <- setNames(color_palette, as.character(seq_len(num_clusters_total)))
metadata(spe)$color_vectors$som_clusters <- cluster_colors

# Genera dinamicamente una palette di colori per i pazienti
num_patients <- length(unique(spe$patient_id))
patient_palette <- scales::hue_pal()(num_patients)
patient_colors <- setNames(patient_palette, as.character(unique(spe$patient_id)))
metadata(spe)$color_vectors$patient_id <- patient_colors

# Costruisci il grafico spaziale usando i cluster Louvain (o manuali)
spe <- buildSpatialGraph(spe, img_id = "patient_id", type = "knn", k = 20)

# Loop attraverso ogni paziente per visualizzare i cluster SOM utilizzando i colori Louvain
for (r in as.vector(unique(spe$patient_id))) {
  p <- plotSpatial(spe[, spe$patient_id == r], 
                   node_color_by = "som_clusters", 
                   img_id = "patient_id", 
                   draw_edges = FALSE, 
                   colPairName = "neighborhood", 
                   nodes_first = FALSE, 
                   node_size_fix = 2,
                   edge_color_fix = "grey") +
    scale_color_manual(values = cluster_colors) +
    ggtitle(paste0("Steinbock interaction graph (SOM): ", r)) + 
    theme(plot.title = element_text(size = 15))
  
  # Salva il grafico per ogni paziente
  ggsave(filename = paste0("Clusters_som_", r, ".png"), p, scale = 1, width = 10, height = 6)
}
