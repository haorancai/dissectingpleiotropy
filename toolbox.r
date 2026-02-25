library(tidygraph)
library(igraph)
# library(MDBED)
library(dplyr)
library(purrr)
bootstrap <- function(tb, n) {
  c(1:n) %>%
    as_tibble() %>%
    mutate(a = map(value, ~ tb %>%
      sample_n(nrow(tb), replace = TRUE) %>%
      mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
      mutate(H = 1.5 * IQR(cor_element)) %>%
      mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
      mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)) %>%
      filter(outlier == FALSE))) %>%
    mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1])) %>%
    unnest(cor_a) %>%
    select(cor_a)
}
bootstrap_with_outlier <- function(tb, n) {
  c(1:n) %>%
    as_tibble() %>%
    mutate(a = map(value, ~ tb %>%
      sample_n(nrow(tb), replace = TRUE) %>%
      mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
      mutate(H = 1.5 * IQR(cor_element)) %>%
      mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
      mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)))) %>%
    mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1])) %>%
    unnest(cor_a) %>%
    select(cor_a)
}


p <- c(0.025, 0.5, 0.975)
p_names <- map_chr(p, ~ paste0(.x * 100, "%"))
p_funs <- map(p, ~ purrr::partial(quantile, probs = .x, na.rm = TRUE)) %>%
  set_names(nm = p_names)


theme_Publication <- function(base_size = 14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size)
  + theme(
      plot.title = element_text(
        face = "bold",
        size = rel(1.2), hjust = 0.5
      ),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = unit(0.2, "cm"),
      legend.margin = unit(0, "cm"),
      legend.title = element_text(face = "italic"),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    ))
}

scale_fill_Publication <- function(...) {
  library(scales)
  discrete_scale("fill", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ...)
}

scale_colour_Publication <- function(...) {
  library(scales)
  discrete_scale("colour", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ...)
}


simulation_for_real_genotype <- function(NrSNP = ncol(genotype), rau = coef,
                                         SNPinfo = genotype, magnitude = 5, mode = "major",
                                         effect_size = FALSE, LD_block, gamma) {
  if (is.null(rau) == TRUE) {
    rau <- runif(1, -1, 1)
  }
  # Specify sample size
  mu <- c(0, 0) # Specify the means of the variables
  sigma <- matrix(c(1, rau, rau, 1), # Specify the covariance matrix of the variables
    ncol = 2
  )

  # effect <- mvrnorm(n = NrSNP, mu = mu, Sigma = sigma)
  effect <- draw.multivariate.laplace(NrSNP, 2, gamma, mu, sigma)
  effect <- as_tibble(cbind(effect, name = colnames(genotype)))


  dummy <- as_tbl_graph(graph_from_data_frame(LD_block %>% dplyr::select(x, y), directed = FALSE, vertices = colnames(X))) %>%
    activate(nodes) %>%
    left_join(effect, by = "name") %>%
    mutate(group = group_components()) %>%
    as_tibble() %>%
    group_by(group) %>%
    nest() %>%
    mutate(data = map(data, function(data) {
      if (nrow(data) > 1) {
        data[1:nrow(data), ] <- data[sample(nrow(data), 1), ]
      }
      data
    })) %>%
    unnest(data) %>%
    ungroup() # %>%
  # dplyr::select(V1,V2)

  # LD_correct_locus <- dummy$name %>% unique()
  effect <- dummy %>%
    dplyr::select(V1, V2) %>%
    mutate_at(c("V1", "V2"), as.numeric)

  effect_size_LD_correct <- dummy %>%
    dplyr::select(-group) %>%
    unique() %>%
    dplyr::select(-name) %>%
    mutate_at(c("V1", "V2"), as.numeric)
  # print(effect_size_LD_correct)
  effect <- as.matrix(effect)


  if (mode == "polygenic") {
    developmental_cor <- cor(effect_size_LD_correct)[2, 1]
    genetic_cor <- cor(SNPinfo %*% effect)[2, 1] # cor(SNPinfo %*% effect)[2,1] / cor(effect_size_LD_correct)[2,1]
    # print(developmental_cor)
    # print(genetic_cor)
  }


  if (mode == "major") {
    big_snp <- round(runif(1, min = 1, max = NrSNP))
    effect[big_snp, ] <- magnitude * effect[big_snp, ]
    developmental_cor <- cor(effect_size_LD_correct)[2, 1]
    genetic_cor <- cor(SNPinfo %*% effect)[2, 1]
  }

  if (mode == "concordant") {
    print("1")
    big_snp <- round(runif(1, min = 1, max = NrSNP))
    effect[big_snp, ] <- magnitude * effect[big_snp, ]


    if (rau > 0) {
      max <- max(abs(effect[big_snp, 1]), abs(effect[big_snp, 2]))
      effect[big_snp, 1] <- max
      effect[big_snp, 2] <- max
    } else {
      max <- max(abs(effect[big_snp, 1]), abs(effect[big_snp, 2]))
      effect[big_snp, 1] <- -max
      effect[big_snp, 2] <- max
    }
    developmental_cor <- cor(effect_size_LD_correct)[2, 1]
    genetic_cor <- cor(SNPinfo %*% effect)[2, 1]
  }


  if (effect_size == TRUE) {
    return(list(effect = effect, genetic_cor = genetic_cor, developmental_cor = developmental_cor))
  } else {
    return(c(genetic_cor, developmental_cor))
  }
}


simulation_for_a_given_population <- function(NrSNP = 1000, rau = coef,
                                              SNPinfo = SNPinfo, magnitude = 30, mode = "major",
                                              effect_size = TRUE) {
  if (is.null(rau) == TRUE) {
    rau <- runif(1, -1, 1)
  }
  # Specify sample size
  mu <- c(0, 0) # Specify the means of the variables
  sigma <- matrix(c(1, rau, rau, 1), # Specify the covariance matrix of the variables
    ncol = 2
  )

  effect <- mvrnorm(n = NrSNP, mu = mu, Sigma = sigma)

  # cor(independent_effect)
  # effect <- rbind(share_effect,independent_effect)

  if (mode == "concordant") {
    big_snp <- round(runif(1, min = 1, max = NrSNP))
    effect[big_snp, ] <- magnitude * effect[big_snp, ]


    if (rau > 0) {
      max <- max(abs(effect[big_snp, 1]), abs(effect[big_snp, 2]))
      effect[big_snp, 1] <- max
      effect[big_snp, 2] <- max
    } else {
      max <- max(abs(effect[big_snp, 1]), abs(effect[big_snp, 2]))
      effect[big_snp, 1] <- -max
      effect[big_snp, 2] <- max
    }
    developmental_cor <- cor(effect[-big_snp, ])[2, 1]
    genetic_cor <- cor(SNPinfo %*% effect)[2, 1]
  }

  if (mode == "polygenic") {
    developmental_cor <- cor(effect)[2, 1]
    genetic_cor <- cor(SNPinfo %*% effect)[2, 1]
  }

  if (mode == "major") {
    big_snp <- round(runif(1, min = 1, max = NrSNP))
    effect[big_snp, ] <- magnitude * effect[big_snp, ]
    developmental_cor <- cor(effect[-big_snp, ])[2, 1]
    genetic_cor <- cor(SNPinfo %*% effect)[2, 1]
  }


  if (effect_size == TRUE) {
    return(list(effect = effect, genetic_cor = genetic_cor, developmental_cor = developmental_cor))
  } else {
    return(c(genetic_cor, developmental_cor))
  }
}




library(LaplacesDemon)
library(MultiRNG)
library(VGAM)

simulation_for_a_given_population_exponential <- function(NrSNP = 300, rau = coef,
                                                          SNPinfo = SNPinfo, effect_size = FALSE, mode, gamma) {
  if (is.null(rau) == TRUE) {
    rau <- runif(1, -1, 1)
  }
  mu <- c(0, 0) # Specify the means of the variables
  sigma <- matrix(c(1, rau, rau, 1), # Specify the covariance matrix of the variables
    ncol = 2
  )

  effect <- draw.multivariate.laplace(NrSNP, 2, gamma, mu, sigma)
  # cor(draw.multivariate.laplace(NrSNP,2, gamma, mu, sigma))
  # tibble(gamma = rgamma(300, 1, scale = 1),laplace =rexp(300, rate = 1)) %>%
  #   #mutate(row = row_number()) %>%
  #   pivot_longer(c(gamma, laplace), names_to = 'type', values_to = 'value') %>%
  #   ggplot(aes(value, fill = type))+geom_density(alpha = .3)
  #
  #

  # if (rau < 0) {
  #   effect[,2] <- -effect[,2]
  # }
  n <- sample(0:10, 1)
  if (n == 0) {
    mode <- "polygenic"
  } else {
    big_snp <- abs(effect) %>%
      as_tibble() %>%
      mutate(num = row_number()) %>%
      mutate(Value = if_else(V1 > V2, V1, V2)) %>%
      arrange(-Value) %>%
      # slice(40) %>%
      top_n(n, Value) %>%
      pull(num)
    # big_snp <- which.max(abs(effect[,2]))
    effect <- as.matrix(effect)
    # print(big_snp)
  }


  if (mode == "concordant") {
    if (rau > 0) {
      max <- apply(rbind(abs(effect[big_snp, 1]), abs(effect[big_snp, 2])), 2, max)
      effect[big_snp, 1] <- max
      effect[big_snp, 2] <- max
    } else {
      max <- apply(rbind(abs(effect[big_snp, 1]), abs(effect[big_snp, 2])), 2, max)
      effect[big_snp, 1] <- -max
      effect[big_snp, 2] <- max
    }

    developmental_cor <- cor(effect[-big_snp, ])[2, 1]
    genetic_cor <- cor(SNPinfo %*% effect)[2, 1]
  }

  if (mode == "polygenic") {
    developmental_cor <- cor(effect)[2, 1]
    genetic_cor <- cor(SNPinfo %*% effect)[2, 1]
  }


  if (mode == "major") {
    max <- apply(rbind(abs(effect[big_snp, 1]), abs(effect[big_snp, 2])), 2, max)
    # Randomly assign each HP locus as concordant (1) or discordant (-1)
    hp_type <- sample(c(1, -1), length(big_snp), replace = TRUE)
    effect[big_snp, 1] <- max
    if (rau > 0) {
      # concordant (hp_type=1): +max,+max (reinforces positive rau)
      # discordant (hp_type=-1): +max,-max (opposes positive rau)
      effect[big_snp, 2] <- max * hp_type
    } else {
      # concordant (hp_type=1): +max,-max (reinforces negative rau)
      # discordant (hp_type=-1): +max,+max (opposes negative rau)
      effect[big_snp, 2] <- -max * hp_type
    }

    developmental_cor <- cor(effect[-big_snp, ])[2, 1]
    genetic_cor <- cor(SNPinfo %*% effect)[2, 1]
  }


  #
  # replace <- sample(NrSNP,20)
  # if (rau > 0){ effect[replace,2] <- effect[replace,1] } else {
  #   effect[replace,2] <- -effect[replace,1]
  # }
  #


  # print('yes')

  if (effect_size == TRUE) {
    return(list(effect = effect, genetic_cor = genetic_cor, developmental_cor = developmental_cor))
  } else {
    return(c(genetic_cor, developmental_cor))
  }
}


get_LOD <- function(cross) {
  # class(cross)[1]<-'dh'
  phen.number <- summary(cross)[[3]]
  # Estimate recombination fractions between markers
  cross <- est.rf(cross)


  newmap <- est.map(cross, error.prob = 0.01)
  cross <- calc.errorlod(cross, error.prob = 0.01)
  top.errorlod(cross)
  cross <- calc.genoprob(cross, step = 1, error.prob = 0.01)
  out.hk <- scanone(cross, pheno.col = 1:phen.number, method = "hk")

  out.hk <- out.hk %>%
    as_tibble() %>%
    mutate(name = rownames(out.hk)) %>%
    dplyr::select(-chr, -pos)
  return(out.hk)
}


get_additive_effect <- function(cross) {
  phen.number <- summary(cross)[[3]]
  cross <- sim.geno(cross, step = 0, n.draws = 16, error.prob = 0.01)

  additive_effect <- c()
  for (i in c(1:phen.number)) {
    a <- effectscan(cross, pheno.col = i, draw = FALSE)
    additive_effect <- bind_cols(additive_effect, a[, 3])
  }

  colnames(additive_effect) <- colnames(cross$pheno)

  return(additive_effect)
}


get_stable_qtl_correlation <- function(cross, k) {
  phen.number <- summary(cross)[[3]]
  cor_additive_effect <- matrix(0, nrow = phen.number, ncol = phen.number)
  cross <- calc.genoprob(cross, step = 1, error.prob = 0.001)
  out.hk <- scanone(cross, pheno.col = 1:phen.number, method = "hk")
  for (i in c(1:k)) {
    cross <- sim.geno(cross, step = 0, n.draws = 16, error.prob = 0.0001)


    additive_effect <- c()
    for (i in c(1:phen.number)) {
      a <- effectscan(cross, pheno.col = i, draw = FALSE)
      additive_effect <- bind_cols(additive_effect, a[, 3])
    }
    cor_additive_effect <- cor_additive_effect + cor(additive_effect)
  }

  colnames(cor_additive_effect) <- colnames(out.hk %>% select(-chr, -pos))
  cor_additive_effect <- cor_additive_effect / k

  return(cor_additive_effect)
}


# Identifing unique loci
findunique <- function(chrlist, collapse = 30) {
  # chrlist is a dataframe with loci position ($pos) and placeholder for unique identifier ($unique, all should be 0), for one chromosome
  while (sum(chrlist$unique == 0) != 0) {
    temp <- which(abs((chrlist$pos - chrlist$pos[which(chrlist$unique == 0)[1]])) <= collapse)
    if (sum(chrlist$unique[temp] != 0) > 0) {
      print("warning")
    }
    # warning means a loci is within collapse distence of two differnt loci (which are not within collapse distence) - look at positions and determine best collapse on your own
    chrlist$unique[temp] <- max(chrlist$unique) + 1
  }
  return(chrlist)
}

# Calling QTL
getQTL <- function(lodcolumn, out, operm, out2, operm2, alpha, collapse = 30) {
  # out,operm,out2 and operm2 are result of scanone and scantwo, lodcolumn is number, alpha is a vector of 2 (for additive and interactions) #example:lodcolumn=1;out=out.hk.L;operm=operm.hk.L;out2=out2.hk.L;operm2=operm2.hk.L;alpha=c(0.1,0.1);cross=cross.L;collapse=30
  # single QTL
  temp1 <- summary(out, perms = operm, alpha = alpha[1], lodcolumn = lodcolumn)
  # addtional QTL on same chromosome from 2d scan
  temp2 <- summary(out2, perms = operm2, alpha = c(0, 0, 0, alpha[1], alpha[1]), lodcolumn = lodcolumn, what = "add", allpairs = F)
  # organize - get chr and pos
  data <- data.frame(chr = c(temp1$chr, temp2$chr1, temp2$chr2), pos = c(temp1$pos, temp2$pos1, temp2$pos2))

  # organize - collapse loci (loci found from 2d with loci found from 1d, take mean pos)
  if (nrow(data) > 0) {
    data <- ldply(split(cbind(data, unique = 0), data$chr), findunique, collapse = collapse)
    data <- aggregate(pos ~ chr + unique, data, mean)
    data <- cbind(data[, c("chr", "pos")], name = paste("Q", 1:nrow(data), sep = ""))
    data$name <- as.character(data$name)
  }

  intbound <- nrow(data)
  # interacting QTL
  temp3 <- summary(out2, perms = operm2, alpha = c(1, 0, alpha[2], 0, 0), lodcolumn = lodcolumn, what = "int")
  int <- matrix(, nrow(temp3), 3)
  if (nrow(temp3) > 0) {
    for (i in 1:nrow(temp3)) {
      # is 1st interacting loci already a QTL
      q1 <- which(data$chr == temp3$chr1[i] & data$pos < (temp3$pos1[i] + collapse) & data$pos > (temp3$pos1[i] - collapse))
      if (length(q1) == 0) {
        # 1st interacting loci is new
        data <- rbind(data, data.frame(chr = temp3$chr1[i], pos = temp3$pos1[i], name = paste("QI", nrow(data) + 1, sep = "")))
        q1 <- nrow(data)
      }
      # is 2st interacting loci already a QTL
      q2 <- which(data$chr == temp3$chr2[i] & data$pos < (temp3$pos2[i] + collapse) & data$pos > (temp3$pos2[i] - collapse))
      if (length(q2) == 0) {
        # 2st interacting loci is new
        data <- rbind(data, data.frame(chr = temp3$chr2[i], pos = temp3$pos2[i], name = paste("QI", nrow(data) + 1, sep = "")))
        q2 <- nrow(data)
      }
      # if two interactions share loci, take mean pos
      if (q1 > intbound) {
        data$pos[q1] <- mean(c(data$pos[q1], temp3$pos1[i]))
      }
      if (q2 > intbound) {
        data$pos[q2] <- mean(c(data$pos[q2], temp3$pos2[i]))
      }
      # info on which loci interact
      int[i, ] <- c(q1, q2, temp3$lod.int[i])
    }
  }

  attr(data, "trait") <- colnames(out)[lodcolumn + 2]
  return(list(data, int))
}


get_replicate_additive_effect <- function(cross, k) {
  replicate_additive_effect <- list()
  for (i in c(1:k)) {
    replicate_additive_effect[i] <- list(get_additive_effect(cross))
  }
  return(replicate_additive_effect)
}


get_element_with_min_absolute_value <- function(var) {
  return(var[which(abs(var) == min(abs(var)))[1]])
}


get_abs_min_from_replicate_additive_effect <- function(cor_additive_effect) {
  return(cor_additive_effect %>%
    map(reshape2::melt) %>%
    bind_rows(.id = "id") %>%
    as_tibble() %>%
    group_by(Var1, Var2) %>%
    summarise(value = get_element_with_min_absolute_value(value)) %>%
    reshape2::acast(Var1 ~ Var2))
}


getQTL <- function(lodcolumn, out, operm, out2, operm2, alpha, collapse = 30) {
  # out,operm,out2 and operm2 are result of scanone and scantwo, lodcolumn is number, alpha is a vector of 2 (for additive and interactions) #example:lodcolumn=1;out=out.hk.L;operm=operm.hk.L;out2=out2.hk.L;operm2=operm2.hk.L;alpha=c(0.1,0.1);cross=cross.L;collapse=30
  # single QTL
  temp1 <- summary(out, perms = operm, alpha = alpha[1], lodcolumn = lodcolumn)
  # addtional QTL on same chromosome from 2d scan
  temp2 <- summary(out2, perms = operm2, alpha = c(0, 0, 0, alpha[1], alpha[1]), lodcolumn = lodcolumn, what = "add", allpairs = F)
  # organize - get chr and pos
  data <- data.frame(chr = c(temp1$chr, temp2$chr1, temp2$chr2), pos = c(temp1$pos, temp2$pos1, temp2$pos2))

  # organize - collapse loci (loci found from 2d with loci found from 1d, take mean pos)
  if (nrow(data) > 0) {
    data <- ldply(split(cbind(data, unique = 0), data$chr), findunique, collapse = collapse)
    data <- aggregate(pos ~ chr + unique, data, mean)
    data <- cbind(data[, c("chr", "pos")], name = paste("Q", 1:nrow(data), sep = ""))
    data$name <- as.character(data$name)
  }

  intbound <- nrow(data)
  # interacting QTL
  temp3 <- summary(out2, perms = operm2, alpha = c(1, 0, alpha[2], 0, 0), lodcolumn = lodcolumn, what = "int")
  int <- matrix(, nrow(temp3), 3)
  if (nrow(temp3) > 0) {
    for (i in 1:nrow(temp3)) {
      # is 1st interacting loci already a QTL
      q1 <- which(data$chr == temp3$chr1[i] & data$pos < (temp3$pos1[i] + collapse) & data$pos > (temp3$pos1[i] - collapse))
      if (length(q1) == 0) {
        # 1st interacting loci is new
        data <- rbind(data, data.frame(chr = temp3$chr1[i], pos = temp3$pos1[i], name = paste("QI", nrow(data) + 1, sep = "")))
        q1 <- nrow(data)
      }
      # is 2st interacting loci already a QTL
      q2 <- which(data$chr == temp3$chr2[i] & data$pos < (temp3$pos2[i] + collapse) & data$pos > (temp3$pos2[i] - collapse))
      if (length(q2) == 0) {
        # 2st interacting loci is new
        data <- rbind(data, data.frame(chr = temp3$chr2[i], pos = temp3$pos2[i], name = paste("QI", nrow(data) + 1, sep = "")))
        q2 <- nrow(data)
      }
      # if two interactions share loci, take mean pos
      if (q1 > intbound) {
        data$pos[q1] <- mean(c(data$pos[q1], temp3$pos1[i]))
      }
      if (q2 > intbound) {
        data$pos[q2] <- mean(c(data$pos[q2], temp3$pos2[i]))
      }
      # info on which loci interact
      int[i, ] <- c(q1, q2, temp3$lod.int[i])
    }
  }

  attr(data, "trait") <- colnames(out)[lodcolumn + 2]
  return(list(data, int))
}

linear <- function(x) {
  x
}

square <- function(x) {
  ifelse(x > 0, x**2, -x**2)
}
#
# get_qtl_correlation <- function(additive_effect, phenotype) {
#   phenotype <- phenotype %>% as_tibble() %>%
#     select(Phen1, Phen2) %>%
#     mutate(cor_a = map2_dbl(Phen1, Phen2, ~ cor(additive_effect %>% select(.x), additive_effect %>% select(.y))[1])) %>%
#     mutate(r_2_a = cor_a ** 2)
#   return(phenotype)}


get_element_squared <- function(m) {
  m <- m %>%
    reshape2::melt() %>%
    mutate(value = ifelse(value > 0, value**2, -(value**2)))


  m <- as.matrix(reshape2::acast(m, m[, 1] ~ m[, 2]))

  return(m)
}


# matrix comparison -------------------------------------------------------

propGmaxAngle <- function(bentE1, bentE2) {
  # Proportion of variation along gmax
  # (Kirkpatrick 2009)
  gmax1 <- bentE1$values[1]
  gmax2 <- bentE2$values[1]
  propGmaxE1 <- gmax1 / sum(bentE1$values)
  propGmaxE2 <- gmax2 / sum(bentE2$values)
  propGmax <- abs(propGmaxE1 - propGmaxE2)
  # propGmaxDiff <- propGmax/((propGmaxE1+propGmaxE2)/2)


  # Angle between the first eigenvector of the two matrices (Ingleby et al. 2014)
  e1gmax <- bentE1$vectors[, 1]
  e2gmax <- bentE2$vectors[, 1]
  angle <- (180 / pi) * acos((e1gmax %*% e2gmax) / (sqrt(e1gmax %*% e1gmax) * sqrt(e2gmax %*% e2gmax)))
  angle <- ifelse(test = angle > 90, yes = 180 - angle, no = angle)

  return(c(propGmax, angle))
}
