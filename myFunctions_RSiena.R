# 1. Functions made by Snijders ----

# To help chose the right n3
source("https://www.stats.ox.ac.uk/~snijders/siena/SE_checks.R")

# To have useful functions for the 5-parameter model for selection (Snijders & Lomi, 2019)
source("https://www.stats.ox.ac.uk/~snijders/siena/SelectionTables.r")

# Corriged selectionTable.norm function
selectionTable.norm <- function(
    x,
    xd,
    name = attr(xd$depvars, "name")[1],
    vname,
    nfirst = x$nwarm + 1,
    multiplier = 1
) {
  
  # Calcul de la base de la sélection
  stab <- selectionTable.basis(
    x, xd, name, vname, 1:2,
    nfirst = nfirst,
    multiplier
  )
  
  vtheta <- stab$vtheta
  vmean <- stab$vmean
  veff <- stab$veff
  veff.eval <- stab$veff.eval   # <-- ligne manquante
  
  # Vérification de egoXaltX
  if (vtheta["egoXaltX"] != 0) {
    cat(
      "Warning: effect of egoXaltX is ",
      vtheta["egoXaltX"],
      ", not equal to 0; \n"
    )
    cat(
      " this function does not apply to such a sienaFit object.\n"
    )
  }
  
  # Gradient
  grad <- matrix(0, length(x$theta), 1)
  
  grad[veff.eval["altX"], 1] <-
    -1 / (2 * vtheta["altSqX"])
  
  grad[veff.eval["altSqX"], 1] <-
    vtheta["altX"] /
    (2 * vtheta["altSqX"] * vtheta["altSqX"])
  
  # Matrice variance-covariance
  covtheta <- x$covtheta
  covtheta[is.na(covtheta)] <- 0
  
  # Vérification de la courbure
  if (vtheta["altSqX"] > 0) {
    cat(
      "Warning: the coefficient of alter squared is positive.\n"
    )
    cat(
      "This means the extremum of the attraction function is a minimum,\n"
    )
    cat(
      "and cannot be interpreted as a social norm.\n"
    )
    flush.console()
  }
  
  list(
    vnorm =
      vmean -
      (vtheta["altX"] /
         (2 * vtheta["altSqX"])),
    
    se.vnorm =
      sqrt(t(grad) %*% covtheta %*% grad)
  )
}

# 2. Functions for GOF plot ----

## 2.1. Geodesic distribution GOF ----

GeodesicDistribution <- function (i, data, sims, period, groupName,
                                  varName, levls = c(1:5, Inf), cumulative = TRUE, ...) {
  # Force chaque processeur à charger RSiena
  require(RSiena) 
  require(sna)    
  
  # Le reste de la fonction ne change pas
  x <- networkExtraction(i, data, sims, period, groupName, varName)
  a <- sna::geodist(symmetrize(x))$gdist
  if (cumulative) {
    gdi <- sapply(levls, function(i){ sum(a <= i) })
  } else {
    gdi <- sapply(levls, function(i){ sum(a == i) })
  }
  names(gdi) <- as.character(levls)
  gdi
}

## 2.2. ESP (OTP) distribution GOF ----

EspDistribution <- function(i, data, sims, period, groupName, varName,
                            levls = 0:20, cumulative = FALSE, ...) {
  require(Matrix)
  
  x <- RSiena::sparseMatrixExtraction(i, data, sims, period, groupName, varName)
  x <- as.matrix(x)
  diag(x) <- 0  # sécurité
  
  # OTP : #{k : i→k et k→j} pour chaque arête i→j
  shared_partners <- (x %*% x)[x == 1]
  
  if (cumulative) {
    esp <- sapply(levls, function(l) sum(shared_partners <= l))  # CDF standard
  } else {
    esp <- sapply(levls, function(l) sum(shared_partners == l))
  }
  
  names(esp) <- as.character(levls)
  return(esp)
}

## 2.3. Automatic GOF function ----

# Pour centrer et scaler (à cause de autograph)
center_scale_gof <- function(gof_obj) {
  gof_scaled <- gof_obj
  for (w in seq_along(gof_obj)) {
    sim <- gof_obj[[w]]$Simulations    # ← RSiena 1.6
    obs <- gof_obj[[w]]$Observations   # ← RSiena 1.6
    
    if (is.null(dim(sim))) sim <- matrix(sim, nrow = 1)
    
    med  <- apply(sim, 2, median, na.rm = TRUE)
    sig  <- apply(sim, 2, sd,     na.rm = TRUE)
    sig[sig == 0] <- 1
    
    gof_scaled[[w]]$Simulations  <- sweep(sweep(sim, 2, med), 2, sig, "/")
    gof_scaled[[w]]$Observations <- (obs - med) / sig
  }
  return(gof_scaled)
}

# Effectuer tous les GOF
gof_all <- function(mod,
                    depvar_name,
                    behvar_name = "mokken_rapFR",
                    levls_in = c(0:8),
                    levls_out = c(0:8),
                    levls_esp = 0:20,
                    levls_geod = c(1:5, Inf),
                    cluster = cl) {
  
  # GOF indegree
  message("Calcul IndegreeDistribution...")
  indegreeGOF <- test_gof(
    object = mod,
    IndegreeDistribution,
    varName = depvar_name,
    join = TRUE,
    levls = levls_in,
    cumulative = TRUE,
    cluster = cluster
  )
  p1 <- plot(indegreeGOF) + theme_light(base_size = 14)
  
  # GOF outdegree
  message("Calcul OutdegreeDistribution...")
  outdegreeGOF <- test_gof(
    object = mod,
    OutdegreeDistribution,
    varName = depvar_name,
    join = TRUE,
    levls = levls_out,
    cumulative = TRUE,
    cluster = cluster
  )
  p2 <- plot(outdegreeGOF) + theme_light(base_size = 14)
  
  # GOF triad census
  message("Calcul TriadCensus...")
  triadCensusGOF <- test_gof(
    mod,
    TriadCensus,
    varName = depvar_name,
    join = TRUE,
    cluster = cluster)
  
  p3_raw  <- plot(center_scale_gof(triadCensusGOF))
  is_text <- sapply(p3_raw$layers, function(l) inherits(l$geom, "GeomText"))
  built   <- ggplot_build(p3_raw)
  text_data        <- built$data[[which(is_text)[1]]]
  text_data$label  <- format(as.integer(round(triadCensusGOF[[1]]$Observations)),
                             big.mark = " ")
  p3 <- p3_raw 
  p3$layers <- p3_raw$layers[!is_text]
  p3 <- p3 + geom_text(data = text_data, aes(x = x, y = y, label = label),
                       color = "red", size = 5.5, vjust = -0.5) +
    labs(y = "Statistic (centered and scaled)") + theme_light(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11))
  
  # GOF ESP-OTP
  message("Calcul EspDistribution...")
  espGOF <- test_gof(
    object = mod,
    EspDistribution,
    varName = depvar_name,
    join = TRUE,
    levls = levls_esp,
    cluster = cluster
  )
  p4 <- plot(espGOF) + theme_light(base_size = 14)
  
  # GOF distances géodésiques
  message("Calcul GeodesicDistribution...")
  geodesicGOF <- test_gof(
    object = mod,
    GeodesicDistribution,
    varName = depvar_name,
    join = TRUE,
    levls = levls_geod,
    cluster = cluster
  )
  p5 <- plot(geodesicGOF) + theme_light(base_size = 14)
  
  # GOF comportement
  message("Calcul BehaviorDistribution...")
  behaviorGOF <- test_gof(
    object = mod,
    BehaviorDistribution,
    varName = behvar_name,
    join = TRUE,
    cluster = cluster
  )
  p6 <- plot(behaviorGOF) + theme_light(base_size = 14)
  
  # GOF ego-alter behavior
  message("Calcul egoAlterCombi...")
  egoAlterCombiGOF <- test_gof(
    object = mod,
    egoAlterCombi,
    varName = c(depvar_name, behvar_name),
    join = TRUE,
    cluster = cluster,
    trafo = function(x) {
      cut(
        x,
        breaks = c(-Inf, 1, 2, Inf),
        labels = FALSE
      )
    }
  )
  
  p7 <- plot(egoAlterCombiGOF) + theme_light(base_size = 14)
  
  net_GOF <- ((p1 + p2) / (p4 + p5) / (p3))
  beh_GOF <- (p6 / p7)
  
  res <- list("indegreeGOF" = p1,
              "outdegreeGOF" = p2,
              "triadCensusGOF" = p3,
              "espGOF" = p4,
              "geodesicGOF" = p5,
              "behaviorGOF" = p6,
              "egoAlterCombiGOF" = p7,
              "net_GOF" = net_GOF,
              "beh_GOF" = beh_GOF)
  
  print(net_GOF)
  print(beh_GOF)
  return(res)
}

# 3. Function for MEMS ----

## 3.1. Fonctions MEMS macro-statistique ----

# Fonction absdiff pour MEMS (macro-statistique = absdiff_tie)
absdiff_mems <- function(network) {
  thesims <- parent.frame()$sim_nets
  j <- parent.frame()$j
  
  # Extraction du comportement simulé
  y <- thesims$sims[[j]][[1]][[2]][[1]] # sims[[j]][[1]][[idx]][[période]] / si période [[1]] -> vague 2 simulée ; [[2]] -> vague 3
  
  # Extraction de la matrice d'adjacence
  mat <- network::as.sociomatrix(network)
  diag(mat) <- 0
  idx <- which(mat != 0, arr.ind = TRUE)
  if (nrow(idx) == 0) return(NA)
  
  mean(abs(y[idx[, 1]] - y[idx[, 2]]), na.rm = TRUE)
}

# Fonction absdiff pour MEMS (macro-statistique = 1 - absdiff_tie / absdiff_all)
absdiff_ratio <- function(network) {
  thesims <- parent.frame()$sim_nets
  j <- parent.frame()$j
  
  y <- thesims$sims[[j]][[1]][[2]][[1]]
  
  mat <- sna::as.sociomatrix(network)
  tie <- which(mat != 0, arr.ind = TRUE)
  absdiff_tie <- mean(abs(y[tie[,1]] - y[tie[,2]]))
  
  D <- abs(outer(y, y, FUN = "-"))
  diag(D) <- NA
  absdiff_all <- mean(D, na.rm = TRUE)
  
  1 - (absdiff_tie / absdiff_all)
}

# 4. Convergence function ----

sienaToConvergence <- function(
    data,
    effects,
    control_algo,
    control_out,
    ans0 = NULL,
    threshold = 0.25,
    tmax_threshold = 0.10,
    thetaBound = 50,
    nbrNodes = 7,
    maxRuns = 10,
    final_estimation = FALSE,
    final_algo = NULL,
    final_out = NULL,
    ...
) {
  
  numr <- 0
  
  if (final_estimation) {
    
    if (is.null(final_algo)) {
      stop(
        "final_estimation = TRUE, mais final_algo n'est pas fourni."
      )
    }
    
    if (is.null(final_out)) {
      stop(
        "final_estimation = TRUE, mais final_out n'est pas fourni."
      )
    }
  }
  
  if (nbrNodes < 1) {
    stop("nbrNodes doit être >= 1.")
  }
  
  if (maxRuns < 1) {
    stop("maxRuns doit être >= 1.")
  }
  
  # Nombre de paramètres estimés
  # (les paramètres de taux ne sont pas dans effects$include)
  p <- sum(effects$include, na.rm = TRUE)
  
  # n2start recommandé par Snijders :
  # 2 * (p + 7) * 2.52^4 / nbrNodes
  n2start_recommended <- ceiling(
    2 * (p + 7) * 2.52^4 / nbrNodes
  )
  
  # Arrondi pratique
  n2start_recommended <- ceiling(
    n2start_recommended / 100
  ) * 100
  
  message(
    "Nombre de paramètres estimés : ", p,
    "\nValeur de départ recommandée pour n2start : ",
    n2start_recommended
  )
  
  
  ans <- siena(
    data = data,
    effects = effects,
    control_algo = control_algo,
    control_out = control_out,
    batch = TRUE,
    silent = FALSE,
    nbrNodes = nbrNodes,
    returnDeps = TRUE,
    prevAns = ans0,
    thetaBound = thetaBound,
    ...
  )
  
  
  
  repeat {
    
    numr <- numr + 1
    
    tm <- as.numeric(ans$tconv.max)
    tx <- as.numeric(ans$tmax)
    
    cat(
      "\nRun", numr,
      "| tconv.max =", round(tm, 4),
      "| tmax =", round(tx, 4),
      "| thetaBound =", thetaBound,
      "\n"
    )
    
    converged <- (
      tm <= threshold &&
        abs(tx) <= tmax_threshold
    )
    
    if (converged) {
      
      message(
        "Convergence atteinte après ",
        numr,
        " estimation(s).",
        "\n  tconv.max = ", round(tm, 4),
        "\n  tmax = ", round(tx, 4)
      )
      
      break
    }
    
    
    
    if (tm > 10) {
      
      stop(
        "Divergence importante : tconv.max > 10."
      )
    }
    
    
    
    if (numr >= maxRuns) {
      
      stop(
        "Convergence non atteinte après ",
        maxRuns,
        " estimation(s).",
        "\n  Dernier tconv.max = ", round(tm, 4),
        "\n  Dernier tmax = ", round(tx, 4)
      )
    }
    
    
    control_algo$nsub <- 1
    
    # On utilise au minimum la valeur recommandée par Snijders,
    # puis on augmente progressivement à chaque relance.
    control_algo$n2start <- max(
      n2start_recommended,
      1000 + numr * 1000
    )
    
    # n3 augmente progressivement.
    # Une valeur d'au moins 5000 est préférable pour la
    # vérification de convergence lorsque le modèle est difficile.
    control_algo$n3 <- max(
      5000,
      250 * p,
      2000 + numr * 1000
    )
    
    # Première continuation : réduire le gain
    if (numr == 1) {
      control_algo$firstg <- 0.01
    }
    
    cat(
      "  -> continuation avec",
      "nsub =", control_algo$nsub,
      ", n2start =", control_algo$n2start,
      ", n3 =", control_algo$n3,
      ", firstg =", control_algo$firstg,
      "\n"
    )
    
    ans <- siena(
      data = data,
      effects = effects,
      control_algo = control_algo,
      control_out = control_out,
      batch = TRUE,
      silent = FALSE,
      nbrNodes = nbrNodes,
      returnDeps = TRUE,
      prevAns = ans,
      thetaBound = thetaBound,
      ...
    )
  }
  
  if (final_estimation) {
    
    message(
      "\nConvergence obtenue : lancement de l'estimation finale."
    )
    
    ans_final <- siena(
      data = data,
      effects = effects,
      control_algo = final_algo,
      control_out = final_out,
      batch = TRUE,
      silent = FALSE,
      nbrNodes = nbrNodes,
      returnDeps = TRUE,
      prevAns = ans,
      thetaBound = thetaBound,
      ...
    )
    
    
    final_tm <- as.numeric(ans_final$tconv.max)
    final_tx <- as.numeric(ans_final$tmax)
    
    cat(
      "\nEstimation finale :",
      "\n  tconv.max =", round(final_tm, 4),
      "\n  tmax =", round(final_tx, 4),
      "\n"
    )
    
    if (
      final_tm > threshold ||
      abs(final_tx) > tmax_threshold
    ) {
      
      warning(
        "L'estimation finale ne satisfait pas complètement ",
        "les critères de convergence.",
        "\n  tconv.max = ", round(final_tm, 4),
        " (seuil = ", threshold, ")",
        "\n  tmax = ", round(final_tx, 4),
        " (seuil = ", tmax_threshold, ")"
      )
    } else {
      
      message(
        "Estimation finale convergée correctement."
      )
    }
    
    return(ans_final)
  }
  
  return(ans)
}

# 5. Jaccard index ----

# Fonction jaccard
jaccard_index <- function(network_w2, network_w4) {
  mat_w2 <- as.sociomatrix(network_w2)
  mat_w4 <- as.sociomatrix(network_w4)
  tab_change <- table(mat_w2, mat_w4)
  res <- tab_change[4] / (tab_change[2] + tab_change[3] + tab_change[4])
  return(res)
}

# 4. Functions for selection plot ---------------------------------------------

selection_plot <- function(
    x,
    xd,
    name,
    vname,
    levls = 0:4,
    levls.alt = seq(0, 4, length.out = 201),
    base_size = 16,
    title = "",
    subtitle = NULL,
    xlab = "Valeur de la variable chez l'alter",
    ylab = "Fonction de sélection",
    legend_title = "Valeur de l'ego",
    save = FALSE,
    filename = NULL,
    width = 10,
    height = 7,
    dpi = 300
) {
  
  # Vérification de ggplot2
  requireNamespace("ggplot2")
  
  # Création du selection plot de Snijders
  p <- selectionTable.plot(
    x = x,
    xd = xd,
    name = name,
    vname = vname,
    levls = levls,
    levls.alt = levls.alt,
    withMax = TRUE,
    base_size = base_size
  )
  
  # Suppression des points intermédiaires de la grille
  # (on conserve les astérisques des maxima)
  p$layers <- p$layers[-1]
  
  # Mise en forme
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = levls,
      labels = levls,
      limits = range(levls.alt)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab,
      colour = legend_title
    ) +
    ggplot2::theme_light(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = base_size + 2
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = base_size - 3
      ),
      axis.title = ggplot2::element_text(
        size = base_size - 1
      ),
      axis.text = ggplot2::element_text(
        size = base_size - 3
      ),
      legend.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = base_size - 4
      ),
      legend.position = "right"
    )
  
  # Sauvegarde optionnelle
  if (save) {
    
    if (is.null(filename)) {
      stop("Tu dois fournir un nom de fichier avec 'filename' si save = TRUE.")
    }
    
    ggplot2::ggsave(
      filename = filename,
      plot = p,
      width = width,
      height = height,
      units = "in",
      dpi = dpi
    )
  }
  
  return(p)
}
