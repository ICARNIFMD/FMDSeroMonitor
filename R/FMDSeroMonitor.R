#' ###############################################################################
# Program for estimation of state-level sero-prevalence rates  #
# Estimation of national level DIVA-positivity rates FMDV infection           #
###############################################################################
#'
#' @param Prevac A n X 3 data frame (n: number of states in the sample) obtained from sero-
#monitoring before going for vaccination (usually collected at zero day of post-vaccination
#                                         monitoring), where row represents the states (with row names as the names of the states), first
#column represents the total sample collected from each state and second column represents the
#number of protected samples, whose Ab-titre is greater than pre-fixed threshold value (e.g. 1.65 in
#                                                                                       log2 scale).
#' @param Postvac A n X 3 data frame (n: number of states in the sample) obtained from sero-
#monitoring after vaccination (usually collected at 28 day of post-vaccination monitoring), where row
#represents the states (with row names as the names of the states), first column represents the total
#sample collected from each state and second column represents the number of protected samples,
#whose Ab-titre is greater than pre-fixed threshold value (e.g. 1.65 in log2 scale) (state/row names in
#                                                                                    both Prevac and Postvac data frame should match).
#' @param Census_Data Animal census N X 1 data frame (state wise bovine population), where rows
#are states and column is the bovine population (cattle + buffalo) (e.g. N = 28 before 2014 &amp; N= 29
#                                                                   after 2014) (row names of the Census_Data should match with those of Prevac and Postvac).
#'
#' @return
#' @export
#'
#' @examples
SeroMonitor <- function (Prevac, Postvac, Census_Data) {
  if(!is.matrix(Prevac) & !is.data.frame(Prevac) & class(Prevac)[1] != "dgCMatrix")
    stop("Wrong input data type of 'Prevac data'")
  if(!is.matrix(Postvac) & !is.data.frame(Postvac) & class(Postvac)[1] != "dgCMatrix")
    stop("Wrong input data type of 'Postvac data'")

  if(sum(is.na(Prevac)) > 0)
    stop("NAs detected in input 'Prevac' data");gc();
  if(sum(is.na(Postvac)) > 0)
    stop("NAs detected in input 'Postvac' data");gc();

  if(sum(Prevac < 0) > 0)
    stop("Negative values detected in input 'Postvac' data");gc();
  if(sum(Postvac < 0) > 0)
    stop("Negative values detected in input 'Postvac' data");gc();

  if(all(Prevac == 0))
    stop("All elements of input 'Prevac' data are zeros");gc();
  if(all(Postvac == 0))
    stop("All elements of input 'Postvac' data are zeros");gc();

  if(any(colSums(Prevac) == 0))
    warning("column size of zeros detected in 'Sero-monitoring pre-vac Data'");gc();
  if(any(colSums(Postvac) == 0))
    warning("column size of zeros detected in 'Sero-monitoring post vac Data'");gc();

  if(!is.matrix(Census_Data) & !is.data.frame(Census_Data) & class(Census_Data)[1] != "dgCMatrix")
    stop("Wrong input data type of 'Census_Data'")
  if(sum(is.na(Census_Data)) > 0)
    warning("NAs detected in input 'Census_Data'");gc();
  if(sum(na.omit(Census_Data) < 0) > 0)
    stop("Negative values detected in input 'Census_Data'");gc();
  if(all(Census_Data == 0))
    stop("All elements of input 'Census_Data' are zeros");gc();

  if (nrow(Census_Data) < nrow(Prevac))
    stop("No. of states in sample must be less than the total states in prevac samples")
  if (nrow(Census_Data) < nrow(Postvac))
    stop("No. of states in sample must be less than the total states in post vac sample")
  if (any(is.na(match(rownames(Prevac), rownames(Census_Data)))) == TRUE)
    warning("State names in input prevac data do not match with names given in server, please write correctly")
  if (any(is.na(match(rownames(Postvac), rownames(Census_Data)))) == TRUE)
    warning("State names in input post vac data do not match with names given in server, please write correctly")

  #Census_Data <- Census_Data[!is.na(match(rownames(SeroSurvData), rownames(Census_Data))), ]
  ###inputs from SeroSurvData
  #N <- 29                 ####total number of states
  ###########Estimation in prevac samples
  N <- nrow(na.omit(Census_Data))
  n1 <- nrow (Prevac)         #######number of states considered in sampling (prevaccination)

  m1i <- as.numeric(Prevac[,1])  #####number random samples from states
  M.all <- as.numeric(Census_Data [, 1])     # state wise animal census data
  M.all <- na.omit (M.all)
  M.total <- sum(M.all)           ####total bovine population in the country
  M1i <- Census_Data[match(row.names(Prevac), row.names(Census_Data)), ]
  prot.pre.s <- as.numeric(Prevac[, 2])
  ##########Estimator
  p1i <- prot.pre.s / m1i                ######state wise proportion of protected animals
  p1i.perc <- round(p1i * 100, 2)
  q1i <- 1 - p1i                    ####rates of failure

  prot.pre <- round(M1i * p1i)        #####state level protected animals at prevaccination

  samp_var.pre <- (p1i * q1i) / (m1i - 1)         #######sample variance
  Est_Var_pre <- (M1i - m1i) / (M1i * (m1i - 1)) * p1i * q1i  ###(Est value) Variance of the estimator

  stand_err.pre <- abs(sqrt(Est_Var_pre))     ##### Standard error of the estimator at prevaccination
  mar_err.pre <- 1.96 * stand_err.pre         #####margin error of the estimator at prevaccination
  up_CI.pre <- (p1i + mar_err.pre) * 100
  low_CI.pre <- (p1i - mar_err.pre) * 100
  CV.pre <- round(((stand_err.pre / p1i) * 100), 2)
  out.prevac <- cbind(p1i.perc, signif(Est_Var_pre, 3), round(stand_err.pre * 100, 2), round(mar_err.pre *100, 2), CV.pre,
                      round(low_CI.pre, 2), round(up_CI.pre, 2), prot.pre)
  #row.names(out1) <-  row.names(SeroSurvData)
  colnam <- c("Protectd(%)", "Variance", "SE(%)", "ME(%)", "CV(%)","95%CI_Low(%)", "95%CI_Upper(%)", "Protected_Tot")
  StateEst_Prevac <- data.frame(out.prevac, row.names = row.names(Prevac))
  rm(out.prevac)
  colnames(StateEst_Prevac) <- colnam
  rm(colnam)
  #print(out1)
  #class(out1) <- "State Level Estimates"

  ####National level estimates
  nat.prot.totpr <- round ((N / n1) * sum (prot.pre))     #####est. value national protected
  nat.prop.pre <- nat.prot.totpr / M.total
  Perc.pre <- round(nat.prop.pre *100, 2)

  bet_var1 <- 1 / (n1 - 1) * sum((prot.pre - sum(prot.pre) / n1) ^ 2)
  within_var1 <- sum ((M1i * (M1i - m1i) / m1i) * samp_var.pre)

  var.nat.tot.pre <- (N * (N - n1) / n1) * bet_var1 + (N / n1) * within_var1
  SE.nat.tot.pre <- round(abs(sqrt(var.nat.tot.pre)), 2)

  var.nat.prop.pre <- var.nat.tot.pre / M.total^2
  SE.nat.prop.pre <- abs(sqrt(var.nat.prop.pre))
  CV.nat.prop.pre <- round(SE.nat.prop.pre / nat.prop.pre * 100, 2)
  ME.nat.prop.pre <- round(1.96 * SE.nat.prop.pre * 100, 2)
  up.CI.nat.prop.pre <- round((nat.prop.pre + 1.96 * SE.nat.prop.pre) * 100 , 2)
  low.CI.nat.prop.pre <- round((nat.prop.pre - 1.96 * SE.nat.prop.pre)* 100 , 2)

  NationalEst.prevac <- c(round(nat.prop.pre, 2), Perc.pre, round(SE.nat.prop.pre * 100, 2), CV.nat.prop.pre, ME.nat.prop.pre,
                          low.CI.nat.prop.pre, up.CI.nat.prop.pre, round(nat.prot.totpr), round(SE.nat.tot.pre, 2))
  names(NationalEst.prevac) <- c("Prot_Prop", "HerdImmunity(%)", "SE(%)", "CV(%)", "ME(%)", "Low95%CI(%)", "Upper95%CI(%)", "PredProtTotal", "SE_PredTotal")


  ###########Estimation in post vaccination samples
  n2 <- nrow (Postvac)         #######number of states considered in sampling (prevaccination)

  m2i <- as.numeric(Postvac[,1])  #####number random samples from states
  M2i <- Census_Data[match(row.names(Postvac), row.names(Census_Data)), ]
  prot.po.s <- as.numeric(Postvac[, 2])
  ##########Estimator
  p2i <- prot.po.s / m2i                ######state wise proportion of protected animals
  p2i.perc <- round(p2i * 100, 2)
  q2i <- 1 - p2i                    ####rates of unprotected animals

  prot.po <- round(M2i * p2i)        #####state level protected animals at post-vaccination

  samp_var.po <- (p2i * q2i) / (m2i - 1)         #######sample variance
  Est_Var_po <- (M2i - m2i) / (M2i * (m2i - 1)) * p2i * q2i  ###(Est value) Variance of the estimator

  stand_err.po <- abs(sqrt(Est_Var_po))     ##### Standard error of the estimator at prevaccination
  mar_err.po <- 1.96 * stand_err.po         #####margin error of the estimator at prevaccination
  up_CI.po <- (p2i + mar_err.po) * 100
  low_CI.po <- (p2i - mar_err.po) * 100
  CV.po <- round(((stand_err.po / p2i) * 100), 2)
  out.postvac <- cbind(p2i.perc, signif(Est_Var_po, 3), round(stand_err.po * 100, 2), round(mar_err.po *100, 2), CV.po,
                       round(low_CI.po, 2), round(up_CI.po, 2), prot.po)
  #row.names(out1) <-  row.names(SeroSurvData)
  colnam <- c("Protectd(%)", "Variance", "SE(%)", "ME(%)", "CV(%)","95%CI_Low(%)", "95%CI_Upper(%)", "Protected_Tot")
  StateEst_Postvac <- data.frame(out.postvac, row.names = row.names(Postvac))
  rm(out.postvac)
  colnames(StateEst_Postvac) <- colnam
  rm(colnam)
  #print(out1)
  #class(out1) <- "State Level Estimates"

  ####National level estimates
  nat.prot.totpo <- round ((N / n2) * sum (prot.po))     #####est. value national protected
  nat.prop.po <- nat.prot.totpo / M.total
  Perc.po <- round(nat.prop.po *100, 2)

  bet_var2 <- 1 / (n2 - 1) * sum((prot.po - sum(prot.po) / n2) ^ 2)
  within_var2 <- sum ((M2i * (M2i - m2i) / m2i) * samp_var.po)

  var.nat.tot.po <- (N * (N - n2) / n2) * bet_var2 + (N / n1) * within_var2
  SE.nat.tot.po <- round(abs(sqrt(var.nat.tot.po)), 2)

  var.nat.prop.po <- var.nat.tot.po / M.total^2
  SE.nat.prop.po <- abs(sqrt(var.nat.prop.po))
  CV.nat.prop.po <- round(SE.nat.prop.po / nat.prop.po * 100, 2)
  ME.nat.prop.po <- round(1.96 * SE.nat.prop.po * 100, 2)
  up.CI.nat.prop.po <- round((nat.prop.po + 1.96 * SE.nat.prop.po) * 100 , 2)
  low.CI.nat.prop.po <- round((nat.prop.po - 1.96 * SE.nat.prop.po)* 100 , 2)

  NationalEst.povac <- c(round(nat.prop.po, 2), Perc.po, round(SE.nat.prop.po * 100, 2), CV.nat.prop.po, ME.nat.prop.po,
                         low.CI.nat.prop.po, up.CI.nat.prop.po, round(nat.prot.totpo), round(SE.nat.tot.po, 2))
  names(NationalEst.povac) <- c("Prot_Prop", "HerdImmunity(%)", "SE(%)", "CV(%)", "ME(%)", "Low95%CI(%)", "Upper95%CI(%)", "PredProtTotal", "SE_PredTotal")

  NatEstimate <- rbind(NationalEst.prevac, NationalEst.povac)
  row.names(NatEstimate) <- c("Prevaccination", "Postvaccination")

  #Impact analysis
  if(nrow(Prevac) > nrow(Postvac)) Prevac = Prevac[match(row.names(Postvac), row.names(Prevac)), ]
  if(nrow(Postvac) > nrow(Prevac)) Postvac = Postvac[match(row.names(Prevac), row.names(Postvac)), ]
  stat <- pval <- vector(mode = "numeric", length = min(nrow(Prevac), nrow(Postvac)))
  for(i in 1:min(nrow(Prevac), nrow(Postvac))){
    #print(i)
    test <- prop.test(x = c(Prevac[i, 2], Postvac [i, 2]), n = c(Prevac[i, 1], Postvac [i, 1]),
                      alternative = "less")
    stat[i] <- test$statistic
    pval[i] <- test$p.value
  }
  stat <- round(stat, 2)
  p1i.1 <- as.numeric(Prevac[, 2])/as.numeric(Prevac[, 1])
  p1i.1 <- round(p1i.1 *100, 2)
  p2i.1 <- as.numeric(Postvac[, 2])/as.numeric(Postvac[, 1])
  p2i.1 <- round(p2i.1 *100, 2)
  adj.pval <- p.adjust(pval, method = "fdr", length(pval))
  logp <- round(-log10(adj.pval), 2)
  impac <- cbind(p1i.1, p2i.1, stat, pval, adj.pval, logp)
  colnames(impac) <- c("Immunity.Pre(%)", "Immunity.Post(%)", "Test-statistic", "p-value", "Adj.p-value", "-log(pval)")
  row.names(impac) <- row.names(Prevac) <- row.names(Postvac)

  result <- list(StateEst.Prevaccination = StateEst_Prevac, StateEst.Postvaccination = StateEst_Postvac, National.Estimate = NatEstimate,
                 Vaccination.impact = impac)
  #save(result, "Result.pdf")
  return(result)
}
#################Server program
Seromonitor_bg <- function (Prevac, Postvac) {
  result <- SeroMonitor (Prevac =Prevac, Postvac = Postvac, Census_Data = census)

  rtffile <- RTF("Ouput.doc")
  addParagraph(rtffile, "Table 1: State level estimation of herd immunity at prevaccination.\n")
  addTable(rtffile, data.frame(result$StateEst.Prevaccination), font.size=9, row.names=T)
  #done(rtffile)
  addParagraph(rtffile, "\n\nTable 2: State level estimation of herd immunity at postvaccination.\n")
  addTable(rtffile, result$StateEst.Postvaccination, font.size=9, row.names=T)
  addParagraph(rtffile, "\n\nTable 3: National level estimation of herd immunity.\n")
  addTable(rtffile, result$National.Estimate, font.size=9, row.names=T)
  addParagraph(rtffile, "\n\nTable 4: Impact of vaccination.\n")
  addTable(rtffile, result$Vaccination.impact, font.size=9, row.names=T)
  done(rtffile)
  return(result)
}

############Visualization of FMD virus sero-prevalence across states in India Map
Seromonitor.map <- function(Prevac, Postvac, shape = coord){

  if (any(is.na(match(rownames(Prevac), shape$ST_NM))) == TRUE)
    warning("State names in input prevac data do not match with names given in server, please write correctly")
  if (any(is.na(match(rownames(Postvac), shape$ST_NM))) == TRUE)
    warning("State names in input post vac data do not match with names given in server, please write correctly")

  #N <- nrow(na.omit(Census_Data))
  n1 <- nrow (Prevac)         #######number of states considered in sampling (prevaccination)
  n2 <-  nrow (Postvac)
  m1i <- as.numeric(Prevac[, 1])  #####number random samples from states during pre-vaccination
  m2i <- as.numeric(Postvac[, 1])  #####number random samples from states during post-vaccination
  prot.pre <- as.numeric(Prevac[, 2])
  prot.po <- as.numeric(Postvac[, 2])
  ##########Estimator
  p1i <- prot.pre / m1i                ######state wise proportion of protected animals
  p2i <- prot.po / m2i                ######state wise proportion of protected animals
  remove(n1, n2, prot.pre, prot.po, m1i, m2i)
  ############Map for prevaccination
  #shape = readShapeSpatial("Admin2.shp")
  id1 <- match(row.names(Prevac), shape$ST_NM)
  id1 <- na.omit(id1)
  id2 <- match(row.names(Postvac), shape$ST_NM)
  id2 <- na.omit(id2)
  prevac.prop <- vector(length = length(shape$ST_NM), mode = "numeric")
  postvac.prop <- vector(length = length(shape$ST_NM), mode = "numeric")
  prevac.prop[id1] <- p1i
  prevac.prop[-id1] <- NA
  prevac.map <- data.frame(id = shape$ST_NM, prevac.prop)
  postvac.prop[id2] <- p2i
  postvac.prop[-id2] <- NA
  postvac.map <- data.frame(id = shape$ST_NM, postvac.prop)
  #remove(DIVA.prop)

  shap.1 <- fortify(shape, region = "ST_NM")
  merge.shap.prevac <- merge(shap.1, prevac.map, by="id", all.x = TRUE)
  merge.shap.postvac <- merge(shap.1, postvac.map, by="id", all.x = TRUE)
  remove(shap.1)
  prevac.plot.ind <- merge.shap.prevac[order(merge.shap.prevac$order), ]
  postvac.plot.ind <- merge.shap.postvac[order(merge.shap.postvac$order), ]

  map.prevac <- ggplot() +
    geom_polygon(data = prevac.plot.ind,
                 aes(x = long, y = lat, group = group, fill = prevac.prop),
                 color = "black", linewidth = 0.25) +
    coord_map() +
    scale_fill_gradient(name = "Proportion", limits=c(0, max(prevac.prop)), low = 'red3', high = 'darkgreen') +
    labs(title = "Estimated herd immunity (Prevaccination)") +
    xlab('Longitude')+
    ylab('Latitude') +
    theme_bw()

  map.postvac <- ggplot() +
    geom_polygon(data = postvac.plot.ind,
                 aes(x = long, y = lat, group = group, fill = postvac.prop),
                 color = "black", linewidth = 0.25) +
    coord_map() +
    scale_fill_gradient(name = "Proportion", limits=c(0, max(postvac.prop)), low = 'red3', high = 'darkgreen') +
    labs(title = "Estimated herd immunity (Postvaccination)") +
    xlab('Longitude')+
    ylab('Latitude') +
    theme_bw()
  map <- marrangeGrob(grobs=list(map.prevac, map.postvac), nrow=1, ncol=2)
  ggsave("Map.pdf", map)
  return(map)
}
