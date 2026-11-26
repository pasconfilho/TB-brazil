# =============================================================================
# Pulmonary tuberculosis notifications, deaths and treatment outcomes in Brazil
# before, during and after the COVID-19 pandemic, 2015-2025
#
# Reproducible analysis script (R >= 4.2).
#
# What this script does, in order:
#   0.  Settings, packages and helper functions
#   1.  Download (optional) and load Sinan tuberculosis microdata
#   2.  Case definitions, selection flow and validation counts
#   3.  Population denominators (IBGE projections, censuses 2010 and 2022)
#   4.  Annual series, age-standardised rates and annual percent change
#   5.  Monthly series, segmented regression and counterfactual (observed/expected)
#   6.  Sensitivity analyses of the counterfactual
#   7.  Rate ratios between periods (sex, age, race/skin colour, region, state)
#   8.  Characteristics of cases, laboratory indicators and HIV by period
#   9.  Treatment outcomes and factors associated with loss to follow-up
#   10. Municipal and spatial analysis (empirical Bayes, Moran's I, LISA, funnel)
#   11. Tuberculosis mortality (Mortality Information System, SIM)
#   12. Figures
#   13. Session information
#
# Inputs (see README.md):
#   data-raw/tb_sinan_BR_2014_2025.csv.gz   Sinan microdata (downloaded in step 1)
#   data/*.csv                               public aggregated files shipped with the repository
# Outputs:
#   output/tables/*.csv, output/figures/*.png, output/results.rds
#
# All analyses use public, anonymised data (Brazilian National Health Council
# Resolution 510/2016). Random seeds are fixed, so simulation-based intervals are
# reproducible up to Monte Carlo error.
# =============================================================================

# ---- 0. Settings, packages and helpers ---------------------------------------

DOWNLOAD_SINAN <- FALSE    # TRUE: download Sinan files from DATASUS with microdatasus
N_SIM          <- 10000    # simulations for prediction intervals
N_PERM         <- 999      # permutations for local Moran's I
SEED           <- 2026

pkgs <- c("data.table", "MASS", "sandwich", "lmtest", "sf", "spdep", "ggplot2", "scales", "R.utils")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install the missing packages first: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table); library(MASS); library(sandwich); library(lmtest)
  library(sf); library(spdep); library(ggplot2)
})
set.seed(SEED)

# Paths are relative to the repository root. Run the script from that folder.
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("data-raw", showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path("output/tables", paste0(name, ".csv")))

YEARS  <- 2015:2025
PERIOD <- list(Pre = 2015:2019, Pandemic = 2020:2021, Post = 2022:2025)
AGES   <- c("0-14", "15-24", "25-34", "35-44", "45-54", "55-64", "65+")

# State codes (IBGE) and regions
UF <- c("11" = "RO", "12" = "AC", "13" = "AM", "14" = "RR", "15" = "PA", "16" = "AP", "17" = "TO",
        "21" = "MA", "22" = "PI", "23" = "CE", "24" = "RN", "25" = "PB", "26" = "PE", "27" = "AL",
        "28" = "SE", "29" = "BA", "31" = "MG", "32" = "ES", "33" = "RJ", "35" = "SP", "41" = "PR",
        "42" = "SC", "43" = "RS", "50" = "MS", "51" = "MT", "52" = "GO", "53" = "DF")
REGION <- c("1" = "North", "2" = "Northeast", "3" = "Southeast", "4" = "South", "5" = "Centre-West")
uf_region <- setNames(REGION[substr(names(UF), 1, 1)], UF)

# Rate ratio with a Poisson 95% confidence interval: (a/pa) / (b/pb)
rate_ratio <- function(a, pa, b, pb) {
  r <- (a / pa) / (b / pb); se <- sqrt(1 / a + 1 / b)
  data.table(RR = r, lo = r * exp(-1.96 * se), hi = r * exp(1.96 * se))
}

# Wilson score interval for a proportion, returned in percent
wilson <- function(k, n, z = 1.96) {
  p <- k / n; den <- 1 + z^2 / n; centre <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  data.table(pct = 100 * p, lo = 100 * (centre - half), hi = 100 * (centre + half))
}

# Annual percent change by Prais-Winsten regression of log10(rate) on year.
# Iterated Prais-Winsten (the first observation is kept and scaled by sqrt(1 - rho^2)).
prais_winsten_apc <- function(years, rates, max_iter = 100, tol = 1e-8) {
  y <- log10(rates); x <- years - years[1]; X <- cbind(1, x); n <- length(y)
  beta <- qr.solve(X, y); rho <- 0
  for (i in seq_len(max_iter)) {
    e <- y - X %*% beta
    rho_new <- sum(e[-1] * e[-n]) / sum(e[-n]^2)
    rho_new <- max(min(rho_new, 0.99), -0.99)
    ys <- c(sqrt(1 - rho_new^2) * y[1], y[-1] - rho_new * y[-n])
    Xs <- rbind(sqrt(1 - rho_new^2) * X[1, ], X[-1, ] - rho_new * X[-n, ])
    fit <- lm.fit(Xs, ys); beta_new <- fit$coefficients
    if (abs(rho_new - rho) < tol && max(abs(beta_new - beta)) < tol) { beta <- beta_new; rho <- rho_new; break }
    beta <- beta_new; rho <- rho_new
  }
  res <- ys - Xs %*% beta; df <- n - 2
  s2 <- sum(res^2) / df; se <- sqrt(diag(s2 * solve(crossprod(Xs))))[2]
  t <- qt(0.975, df); f <- function(b) (10^b - 1) * 100
  lo <- f(beta[2] - t * se); hi <- f(beta[2] + t * se)
  data.table(APC = f(beta[2]), lo = lo, hi = hi, p = 2 * pt(-abs(beta[2] / se), df), rho = rho,
             trend = fifelse(lo > 0, "Increasing", fifelse(hi < 0, "Decreasing", "Stationary")))
}

# ---- 1. Sinan microdata ------------------------------------------------------

SINAN_FILE <- "data-raw/tb_sinan_BR_2014_2025.csv.gz"
SINAN_VARS <- c("DT_NOTIFIC", "NU_ANO", "SG_UF_NOT", "ID_MUNICIP", "DT_DIAG", "NU_IDADE_N", "CS_SEXO",
                "CS_GESTANT", "CS_RACA", "CS_ESCOL_N", "SG_UF", "ID_MN_RESI", "NDUPLIC_N", "IN_VINCULA",
                "TRATAMENTO", "FORMA", "AGRAVAIDS", "AGRAVALCOO", "AGRAVDIABE", "AGRAVDOENC", "AGRAVDROGA",
                "AGRAVTABAC", "BACILOSC_E", "CULTURA_ES", "TEST_MOLEC", "TEST_SENSI", "HISTOPATOL", "HIV",
                "ANT_RETRO", "TRAT_SUPER", "TRATSUP_AT", "DT_INIC_TR", "SITUA_ENCE", "DT_ENCERRA",
                "POP_LIBER", "POP_RUA", "POP_SAUDE", "POP_IMIG", "BENEF_GOV", "TRANSF", "UF_TRANSF", "MUN_TRANSF")

if (DOWNLOAD_SINAN) {
  # microdatasus >= 3.0.0.9000 (development branch) supports Sinan tuberculosis.
  # Install with: remotes::install_github("rfsaldanha/microdatasus", ref = "dev")
  sinan <- microdatasus::fetch_datasus(year_start = 2014, year_end = 2025,
                                       information_system = "SINAN-TUBERCULOSE",
                                       vars = SINAN_VARS, track_source = TRUE, provenance = TRUE, timeout = 1800)
  fwrite(as.data.table(microdatasus::datasus_provenance(sinan)), "data-raw/provenance_sinan.csv")
  fwrite(sinan, SINAN_FILE)
}
if (!file.exists(SINAN_FILE)) stop("Sinan file not found: set DOWNLOAD_SINAN <- TRUE or place it in data-raw/.")
d <- fread(SINAN_FILE, colClasses = "character", na.strings = c("", "NA"))

# ---- 2. Case definitions and selection flow ----------------------------------

d[, dt_diag := as.IDate(DT_DIAG)]
d[, `:=`(year = year(dt_diag), month = month(dt_diag))]
d[, cuf := substr(ID_MN_RESI, 1, 2)]
d[, `:=`(uf = UF[cuf], region = REGION[substr(cuf, 1, 1)], mun6 = substr(ID_MN_RESI, 1, 6))]

# Laboratory confirmation: positive smear, positive culture, or M. tuberculosis detected
# by the rapid molecular test (TEST_MOLEC 1 = detected rifampicin-sensitive, 2 = rifampicin-resistant)
d[, lab := BACILOSC_E %chin% "1" | CULTURA_ES %chin% "1" | TEST_MOLEC %chin% c("1", "2")]

# Age: NU_IDADE_N uses a unit digit (4 = years; 1-3 = hours/days/months, i.e. < 1 year)
d[, age := {
  u <- substr(NU_IDADE_N, 1, 1); v <- suppressWarnings(as.integer(substr(NU_IDADE_N, 2, 4)))
  fifelse(nchar(NU_IDADE_N) == 4 & u == "4", v, fifelse(nchar(NU_IDADE_N) == 4 & u %chin% c("1", "2", "3"), 0L, NA_integer_))
}]
d[, age_group := as.character(cut(age, c(-1, 14, 24, 34, 44, 54, 64, 200), labels = AGES))]
d[, sex := fcase(CS_SEXO == "M", "Male", CS_SEXO == "F", "Female")]
d[, race := fcase(CS_RACA == "1", "White", CS_RACA == "2", "Black", CS_RACA == "3", "Asian",
                  CS_RACA == "4", "Brown (pardo)", CS_RACA == "5", "Indigenous")]
d[, outcome := suppressWarnings(as.integer(SITUA_ENCE))]
yesno <- function(x) fcase(x == "1", 1L, x == "2", 0L)          # 9 / blank -> NA
d[, `:=`(pdl = yesno(POP_LIBER), homeless = yesno(POP_RUA), alcohol = yesno(AGRAVALCOO),
         drugs = yesno(AGRAVDROGA), tobacco = yesno(AGRAVTABAC), diabetes = yesno(AGRAVDIABE))]
d[, hiv := fcase(HIV == "1", "Positive", HIV == "2", "Negative", HIV == "3", "Pending", HIV == "4", "Not done")]

flow <- data.table(step = "TB records in national Sinan files 2014-2025", n = nrow(d))
s <- d[year %between% c(2015, 2025)];       flow <- rbind(flow, list("Diagnosed 2015-2025", nrow(s)))
s <- s[!is.na(uf)];                          flow <- rbind(flow, list("Resident in Brazil (valid municipality)", nrow(s)))
s <- s[FORMA == "1"];                        flow <- rbind(flow, list("Pulmonary only", nrow(s)));  all_pulm <- s
s <- s[lab == TRUE];                         flow <- rbind(flow, list("Laboratory-confirmed", nrow(s))); lab_pulm <- s
cn <- s[TRATAMENTO == "1"];                  flow <- rbind(flow, list("New cases (study population)", nrow(cn)))
cn[, period := fifelse(year <= 2019, "Pre", fifelse(year <= 2021, "Pandemic", "Post"))]
save_table(flow, "flow")

# Validation against official tabulations (TabNet and the 2026 national bulletin)
validation <- data.table(
  check = c("Laboratory-confirmed pulmonary cases, 2015-2025, all entry types (compare with TabNet)",
            "New cases of all forms in 2025 (compare with the national bulletin, 84,368)"),
  n = c(d[year %between% c(2015, 2025) & FORMA == "1" & lab == TRUE, .N], d[year == 2025 & TRATAMENTO == "1", .N]))
save_table(validation, "validation")
print(flow); print(validation)

# Completeness of the deprivation-of-liberty field and Sinan record-management flags among study
# cases. The public files carry no personal identifiers, so no further record linkage is possible;
# records are counted as released, which is also how the official TabNet tabulations count them.
data_quality <- cn[, .(n = .N, pct_pdl_field_missing = mean(is.na(pdl)) * 100,
                       pct_linked_IN_VINCULA = mean(IN_VINCULA %chin% "1") * 100,
                       pct_NDUPLIC_N_1 = mean(NDUPLIC_N %chin% "1") * 100), keyby = year]
save_table(data_quality, "data_quality_by_year")

# Entry type (TRATAMENTO) among all laboratory-confirmed pulmonary cases, by year. Only code 1
# (new case) enters the study population; the other codes show where previously treated cases come from.
ENTRY <- c("1" = "New", "2" = "Relapse", "3" = "Re-entry after loss to follow-up", "4" = "Unknown",
           "5" = "Transfer", "6" = "Post-death")
entry_types <- lab_pulm[, .N, keyby = .(year, entry = ENTRY[TRATAMENTO])][, pct := N / sum(N) * 100, by = year]
entry_types <- dcast(entry_types, year ~ entry, value.var = "pct", fill = 0)
save_table(entry_types, "entry_types_by_year")

# ---- 3. Population denominators ----------------------------------------------

# IBGE population projections, 2024 revision: one row per state (or BR) and group
# (7 age groups plus "Homens"/"Mulheres"), columns 2015-2025.
P <- fread("data/pop_uf_ibge_projections_2024.csv", header = TRUE)
popget <- function(u, g) unlist(P[uf == u & grupo == g, as.character(YEARS), with = FALSE])
pop_br     <- popget("BR", "Homens") + popget("BR", "Mulheres")
pop_sex_br <- list(Male = popget("BR", "Homens"), Female = popget("BR", "Mulheres"))
pop_age_br <- sapply(AGES, function(a) popget("BR", a))                       # 11 x 7
pop_uf     <- sapply(UF, function(u) popget(u, "Homens") + popget(u, "Mulheres"))  # 11 x 27
rownames(pop_uf) <- rownames(pop_age_br) <- YEARS; colnames(pop_uf) <- UF
pop_region <- sapply(REGION, function(r) rowSums(pop_uf[, uf_region == r, drop = FALSE]))
colnames(pop_region) <- REGION
std_weights <- pop_age_br["2022", ] / sum(pop_age_br["2022", ])                # 2022 Brazilian population

# Race/skin colour: the state race composition is interpolated linearly between the 2010 and
# 2022 censuses, kept at the 2022 composition for 2023-2025, and applied to annual state
# populations. Self-identification as Black or brown rose between censuses, so a fixed 2022
# denominator would understate past denominators for these groups (kept as a sensitivity analysis).
# Census columns (IBGE categories branca, preta, amarela, parda, indigena) and labels used in the tables
RACES <- c(White = "White", Black = "Black", Asian = "Asian", Brown = "Brown (pardo)", Indigenous = "Indigenous")
r10 <- fread("data/race_uf_census2010.csv", header = TRUE, colClasses = list(character = "cuf"))
r22 <- fread("data/race_uf_census2022.csv", header = TRUE, colClasses = list(character = "cuf"))
share <- function(r) { m <- as.matrix(r[, names(RACES), with = FALSE]); rownames(m) <- r$cuf; m / rowSums(m) }
s10 <- share(r10); s22 <- share(r22)
pop_race <- matrix(0, length(YEARS), 5, dimnames = list(YEARS, RACES))
for (cu in names(UF)) for (y in YEARS) {
  w <- min((y - 2010) / 12, 1)
  pop_race[as.character(y), ] <- pop_race[as.character(y), ] + ((1 - w) * s10[cu, ] + w * s22[cu, ]) * pop_uf[as.character(y), UF[cu]]
}
pop_race_fixed <- colSums(as.matrix(r22[, names(RACES), with = FALSE])); names(pop_race_fixed) <- RACES
# Sensitivity: self-identification may have kept changing after 2022, so the 2010-2022 linear trend in
# each state's composition is also extended to 2023-2025 (weights above 1).
pop_race_extrap <- matrix(0, length(YEARS), 5, dimnames = list(YEARS, RACES))
for (cu in names(UF)) for (y in YEARS) {
  w <- (y - 2010) / 12
  pop_race_extrap[as.character(y), ] <- pop_race_extrap[as.character(y), ] + ((1 - w) * s10[cu, ] + w * s22[cu, ]) * pop_uf[as.character(y), UF[cu]]
}

# ---- 4. Annual series --------------------------------------------------------

# Count of cases per year, as a named vector 2015-2025 (zeros where no cases)
by_year <- function(x) {
  v <- x[, .N, by = year]; out <- setNames(v$N[match(YEARS, v$year)], YEARS); out[is.na(out)] <- 0; out
}
age_std_rate <- function(x, pop_age) {
  cnt <- dcast(x[!is.na(age_group), .N, by = .(year, age_group)], year ~ age_group, value.var = "N", fill = 0)
  m <- as.matrix(cnt[match(YEARS, year), AGES, with = FALSE]); m[is.na(m)] <- 0
  as.vector((m / pop_age * 1e5) %*% std_weights)
}
annual <- data.table(year = YEARS, cases = by_year(cn), pop = pop_br, notif_lab = by_year(lab_pulm),
                     pulm_new_all = by_year(all_pulm[TRATAMENTO == "1"]))
annual[, `:=`(rate = cases / pop * 1e5, asr = age_std_rate(cn, pop_age_br),
              rate_male = by_year(cn[sex == "Male"]) / pop_sex_br$Male * 1e5,
              rate_female = by_year(cn[sex == "Female"]) / pop_sex_br$Female * 1e5,
              pdl = by_year(cn[pdl == 1]))]
annual[, `:=`(pdl_pct = pdl / cases * 100, rate_non_pdl = (cases - pdl) / pop * 1e5,
              pct_not_new = (1 - cases / notif_lab) * 100,          # all entry types other than "new case"
              pct_lab_among_new_pulmonary = by_year(all_pulm[TRATAMENTO == "1" & lab == TRUE]) / pulm_new_all * 100)]
xpert <- cn[, .(pct_xpert_done = mean(TEST_MOLEC %chin% c("1", "2", "3", "4")) * 100,
                pct_xpert_only = mean(TEST_MOLEC %chin% c("1", "2") & !BACILOSC_E %chin% "1" & !CULTURA_ES %chin% "1") * 100), by = year]
annual <- merge(annual, xpert, by = "year")
save_table(annual, "annual_series")

apc_rows <- list()
apc_series <- c("Brazil, crude" = "rate", "Brazil, age-standardised" = "asr", "Brazil, excluding PDL" = "rate_non_pdl",
                "Brazil, male" = "rate_male", "Brazil, female" = "rate_female")
region_rate <- sapply(REGION, function(r) by_year(cn[region == r]) / pop_region[, r] * 1e5)
colnames(region_rate) <- REGION
# Windows: pre-pandemic, post-pandemic (2022-2025, without the pandemic years) and the whole period.
for (nm in names(apc_series)) for (win in list(2015:2019, 2022:2025, 2015:2025))
  apc_rows[[length(apc_rows) + 1]] <- cbind(series = nm, window = paste(range(win), collapse = "-"),
                                            prais_winsten_apc(win, annual[year %in% win][[apc_series[nm]]]))
for (r in REGION) for (win in list(2015:2019, 2022:2025, 2015:2025))
  apc_rows[[length(apc_rows) + 1]] <- cbind(series = r, window = paste(range(win), collapse = "-"),
                                            prais_winsten_apc(win, region_rate[as.character(win), r]))
apc <- rbindlist(apc_rows); save_table(apc, "table_s2_apc")

# ---- 5. Monthly series, segmented regression and counterfactual ---------------

months <- seq(as.IDate("2015-01-01"), as.IDate("2025-12-01"), by = "month")
T_PAND <- 62   # index (0-based) of March 2020
T_POST <- 84   # index (0-based) of January 2022

# Monthly frame: count y, population interpolated between mid-year (July) estimates with
# linear extrapolation at the ends, Fourier terms for seasonality and segment terms.
monthly_frame <- function(counts, pop_annual) {
  m <- data.table(date = months, t = 0:131)
  m[, `:=`(year = year(date), month = month(date))]
  m[, y := as.numeric(counts)]
  idx <- (m$year - 2015) * 12 + m$month - 1; mid <- (YEARS - 2015) * 12 + 6; pv <- as.numeric(pop_annual)
  s0 <- (pv[2] - pv[1]) / 12; s1 <- (pv[11] - pv[10]) / 12
  pm <- approx(mid, pv, xout = idx, rule = 2)$y
  pm[idx < mid[1]] <- pv[1] + (idx[idx < mid[1]] - mid[1]) * s0
  pm[idx > mid[11]] <- pv[11] + (idx[idx > mid[11]] - mid[11]) * s1
  m[, `:=`(pop = pm, days = days_in_month)]
  for (k in 1:2) m[, (paste0(c("s", "c"), k)) := list(sin(2 * pi * k * month / 12), cos(2 * pi * k * month / 12))]
  m[, `:=`(pand = as.integer(t >= T_PAND & t < T_POST), post = as.integer(t >= T_POST),
           t_pand = pmax(t - T_PAND, 0) * (t >= T_PAND), t_post = pmax(t - T_POST, 0) * (t >= T_POST))]
  m[, off := log(pop) + log(days)]
  m[]
}
days_in_month <- as.numeric(diff(c(months, as.IDate("2026-01-01"))))   # 28-31 days
# Count of cases per calendar month, January 2015 to December 2025 (132 values)
monthly_counts <- function(x) {
  v <- x[, .N, by = .(year, month)]
  out <- v$N[match(paste(year(months), month(months)), paste(v$year, v$month))]; out[is.na(out)] <- 0; out
}

# Counterfactual: quasi-Poisson with linear trend and two Fourier pairs fitted to the
# pre-pandemic months (from `start` to February 2020), projected forward. Prediction
# intervals: coefficients drawn from their multivariate normal distribution and counts drawn
# from a gamma-Poisson mixture with the estimated overdispersion.
counterfactual <- function(m, trend = TRUE, start = 0, end = 132, nsim = N_SIM, return_sims = FALSE) {
  terms <- c(if (trend) "t", "s1", "c1", "s2", "c2")
  pre <- m[t >= start & t < T_PAND]
  fit <- glm(reformulate(terms, "y"), family = quasipoisson(), offset = off, data = pre)
  phi <- max(summary(fit)$dispersion, 1.0001)
  X <- cbind(1, as.matrix(m[, ..terms])); b <- coef(fit); V <- vcov(fit)
  expected <- as.vector(exp(X %*% b + m$off))
  B <- mvrnorm(nsim, b, V); MU <- exp(B %*% t(X) + matrix(m$off, nsim, nrow(m), byrow = TRUE))
  YS <- matrix(rpois(length(MU), rgamma(length(MU), shape = MU / (phi - 1), scale = phi - 1)), nsim)
  windows <- list(Pandemic = which(m$t >= T_PAND & m$t < T_POST), Post = which(m$t >= T_POST & m$t < end),
                  Cumulative = which(m$t >= T_PAND & m$t < end))
  for (y in 2020:2025) windows[[as.character(y)]] <- which(m$year == y & m$t >= T_PAND)
  out <- rbindlist(lapply(names(windows), function(w) {
    i <- windows[[w]]; o <- sum(m$y[i]); e <- sum(expected[i]); s <- rowSums(YS[, i, drop = FALSE])
    data.table(window = w, observed = o, expected = e, difference = o - e,
               OE = o / e, OE_lo = o / quantile(s, 0.975), OE_hi = o / quantile(s, 0.025))
  }))
  res <- list(summary = out, expected = expected, pi_lo = apply(YS, 2, quantile, 0.025),
              pi_hi = apply(YS, 2, quantile, 0.975), fit = fit)
  if (return_sims) res$sims <- YS
  res
}

# Records with the deprivation-of-liberty field blank are counted outside prisons (non-PDL).
cases_monthly <- list(all = monthly_counts(cn), non_pdl = monthly_counts(cn[is.na(pdl) | pdl != 1]),
                      pdl = monthly_counts(cn[pdl == 1]))
mf <- lapply(cases_monthly, monthly_frame, pop_annual = pop_br)

# Segmented (interrupted time-series) quasi-Poisson model with Newey-West standard errors
its <- glm(y ~ t + pand + t_pand + post + t_post + s1 + c1 + s2 + c2, family = quasipoisson(),
           offset = off, data = mf$all)
V_nw <- NeweyWest(its, lag = 12, prewhite = FALSE)
its_tab <- as.data.table(coeftest(its, vcov. = V_nw)[, ], keep.rownames = "term")
annual_change <- function(terms) {
  b <- sum(coef(its)[terms]); se <- sqrt(sum(V_nw[terms, terms]))
  f <- function(v) (exp(12 * v) - 1) * 100
  data.table(pct_per_year = f(b), lo = f(b - 1.96 * se), hi = f(b + 1.96 * se))
}
its_trends <- rbind(cbind(segment = "Pre-pandemic", annual_change("t")),
                    cbind(segment = "Pandemic", annual_change(c("t", "t_pand"))),
                    cbind(segment = "Post-pandemic", annual_change(c("t", "t_pand", "t_post"))))
its_levels <- its_tab[term %in% c("pand", "post"), .(term, IRR = exp(Estimate), lo = exp(Estimate - 1.96 * `Std. Error`),
                                                      hi = exp(Estimate + 1.96 * `Std. Error`), p = `Pr(>|z|)`)]
save_table(its_tab, "its_coefficients"); save_table(its_trends, "its_annual_trends"); save_table(its_levels, "its_level_changes")

cf <- lapply(mf, counterfactual)
cf_table <- rbindlist(lapply(names(cf), function(k) cbind(series = k, cf[[k]]$summary)))
save_table(cf_table, "table3_counterfactual_brazil")

# By region and state, excluding PDL (and all cases)
cf_geo <- list()
for (r in REGION) for (who in c("all", "non_pdl")) {
  x <- cn[region == r]; if (who == "non_pdl") x <- x[is.na(pdl) | pdl != 1]
  cf_geo[[length(cf_geo) + 1]] <- cbind(level = "Region", unit = r, series = who,
                                        counterfactual(monthly_frame(monthly_counts(x), pop_region[, r]), nsim = 5000)$summary)
}
for (u in UF) for (who in c("all", "non_pdl")) {
  x <- cn[uf == u]; if (who == "non_pdl") x <- x[is.na(pdl) | pdl != 1]
  cf_geo[[length(cf_geo) + 1]] <- cbind(level = "State", unit = u, series = who,
                                        counterfactual(monthly_frame(monthly_counts(x), pop_uf[, u]), nsim = 5000)$summary)
}
cf_geo <- rbindlist(cf_geo); save_table(cf_geo, "table_s3_counterfactual_regions_states")

# ---- 6. Sensitivity analyses of the counterfactual ---------------------------

# Negative binomial baseline (theta by maximum likelihood) with the same simulation logic
counterfactual_nb <- function(m, nsim = N_SIM) {
  pre <- m[t < T_PAND]
  fit <- glm.nb(y ~ t + s1 + c1 + s2 + c2 + offset(off), data = pre)
  X <- cbind(1, as.matrix(m[, .(t, s1, c1, s2, c2)])); b <- coef(fit); th <- fit$theta
  expected <- as.vector(exp(X %*% b + m$off))
  MU <- exp(mvrnorm(nsim, b, vcov(fit)) %*% t(X) + matrix(m$off, nsim, nrow(m), byrow = TRUE))
  YS <- matrix(rpois(length(MU), rgamma(length(MU), shape = th, scale = MU / th)), nsim)
  rbindlist(lapply(list(Pandemic = which(m$t >= T_PAND & m$t < T_POST), Post = which(m$t >= T_POST),
                        Cumulative = which(m$t >= T_PAND)), function(i) {
    o <- sum(m$y[i]); s <- rowSums(YS[, i, drop = FALSE])
    data.table(observed = o, expected = sum(expected[i]), OE = o / sum(expected[i]),
               OE_lo = o / quantile(s, 0.975), OE_hi = o / quantile(s, 0.025))
  }), idcol = "window")
}

# Moving-block bootstrap: 12-month blocks of the pre-pandemic log residuals are resampled and added
# to the projected log expected counts, so the prediction interval keeps any residual autocorrelation.
block_bootstrap_pi <- function(m, B = 4000, L = 12) {
  res <- counterfactual(m, nsim = 2000); fit <- res$fit
  pre <- m$t < T_PAND
  lr <- log(pmax(m$y[pre], 0.5) / res$expected[pre]); lr <- lr - mean(lr)
  X <- cbind(1, as.matrix(m[, .(t, s1, c1, s2, c2)])); Bc <- mvrnorm(B, coef(fit), vcov(fit))
  nb <- ceiling(nrow(m) / L); post <- m$t >= T_POST
  sp <- vapply(seq_len(B), function(b) {
    st <- sample.int(length(lr) - L + 1, nb, replace = TRUE)
    eps <- unlist(lapply(st, function(s) lr[s:(s + L - 1)]))[seq_len(nrow(m))]
    sum(exp(X[post, ] %*% Bc[b, ] + m$off[post] + eps[post]))
  }, numeric(1))
  r <- residuals(fit, type = "pearson")
  o <- sum(m$y[post])
  data.table(rho1 = cor(r[-1], r[-length(r)]), ljung_box_12_p = Box.test(r, lag = 12, type = "Ljung-Box")$p.value,
             OE_post = res$summary[window == "Post", OE],
             PI_param_lo = res$summary[window == "Post", OE_lo], PI_param_hi = res$summary[window == "Post", OE_hi],
             PI_block_lo = o / quantile(sp, 0.975), PI_block_hi = o / quantile(sp, 0.025))
}

deaths_monthly <- fread("data/sim_tb_deaths_monthly_uf.csv", header = TRUE)[order(year, month)]
mf$deaths <- monthly_frame(deaths_monthly$BR, pop_br)

sens <- list()
for (k in c("all", "non_pdl", "pdl")) {
  m <- mf[[k]]
  specs <- list("Main: quasi-Poisson, trend + seasonality, Jan 2015-Feb 2020" = counterfactual(m)$summary,
                "Shorter baseline (Mar 2017-Feb 2020), trend" = counterfactual(m, start = 26)$summary,
                "No trend, baseline Mar 2017-Feb 2020 (flat rate)" = counterfactual(m, trend = FALSE, start = 26)$summary,
                "No trend, baseline Mar 2018-Feb 2020 (flat rate)" = counterfactual(m, trend = FALSE, start = 38)$summary,
                "Negative binomial, trend + seasonality" = counterfactual_nb(m),
                "Post-pandemic window 2022-2024 (2025 data preliminary)" = counterfactual(m, end = 120)$summary)
  sens[[k]] <- rbindlist(lapply(specs, function(x) x[window %in% c("Pandemic", "Post", "Cumulative"),
                                                     .(window, observed, expected, OE, OE_lo, OE_hi)]), idcol = "specification")
}
# Excluding states with known data-quality problems (non-PDL cases, population without those states)
for (drop in list("PB", c("PB", "PE", "DF"))) {
  x <- cn[!uf %chin% drop & (is.na(pdl) | pdl != 1)]
  m <- monthly_frame(monthly_counts(x), pop_br - rowSums(pop_uf[, drop, drop = FALSE]))
  sens[[paste0("non_pdl_excl_", paste(drop, collapse = "_"))]] <-
    cbind(specification = paste("Excluding", paste(drop, collapse = ", ")),
          counterfactual(m)$summary[window %in% c("Pandemic", "Post", "Cumulative"), .(window, observed, expected, OE, OE_lo, OE_hi)])
}
# Alternative case definitions (all new pulmonary cases; all laboratory-confirmed cases, new and retreatment)
alt_def <- list("All new pulmonary cases (any confirmation)" = all_pulm[TRATAMENTO == "1"],
                "All laboratory-confirmed pulmonary cases (all entry types)" = lab_pulm)
sens$alternative_definitions <- rbindlist(lapply(alt_def, function(x)
  counterfactual(monthly_frame(monthly_counts(x), pop_br))$summary[window %in% c("Pandemic", "Post", "Cumulative"),
                                                                   .(window, observed, expected, OE, OE_lo, OE_hi)]),
  idcol = "specification")
sens_table <- rbindlist(sens, idcol = "series", fill = TRUE)
save_table(sens_table, "table_s1_sensitivity")

autocorr <- rbindlist(lapply(c(all = "all", non_pdl = "non_pdl", pdl = "pdl", deaths = "deaths"),
                             function(k) block_bootstrap_pi(mf[[k]])), idcol = "series")
save_table(autocorr, "table_s1_block_bootstrap")

# ---- 7. Rate ratios between periods ------------------------------------------

rr_row <- function(group, category, counts, pop) {
  o <- data.table(group = group, category = category)
  for (p in names(PERIOD)) {
    yy <- as.character(PERIOD[[p]])
    o[, (paste0("n_", p)) := sum(counts[yy])]
    o[, (paste0("rate_", p)) := sum(counts[yy]) / sum(pop[yy]) * 1e5]
  }
  pre <- as.character(PERIOD$Pre)
  for (p in c("Pandemic", "Post")) {
    yy <- as.character(PERIOD[[p]])
    r <- rate_ratio(sum(counts[yy]), sum(pop[yy]), sum(counts[pre]), sum(pop[pre]))
    o[, (paste0("RR_", p, c("", "_lo", "_hi"))) := r]
  }
  o
}
rr <- list(rr_row("Total", "Brazil", by_year(cn), pop_br),
           rr_row("Total", "Excluding PDL", by_year(cn[is.na(pdl) | pdl != 1]), pop_br))
for (sx in c("Male", "Female")) rr[[length(rr) + 1]] <- rr_row("Sex", sx, by_year(cn[sex == sx]), pop_sex_br[[sx]])
for (a in AGES) rr[[length(rr) + 1]] <- rr_row("Age group (years)", a, by_year(cn[age_group == a]), pop_age_br[, a])
# Asian is reported in no table: the Asian population fell by more than half between censuses
# (2.1 to 0.85 million), which points to reclassification, so its rates are not interpretable.
for (rc in setdiff(RACES, "Asian"))
  rr[[length(rr) + 1]] <- rr_row("Race/skin colour", rc, by_year(cn[race == rc]), pop_race[, rc])
for (r in REGION) rr[[length(rr) + 1]] <- rr_row("Region", r, by_year(cn[region == r]), pop_region[, r])
rr <- rbindlist(rr); save_table(rr, "table2_rate_ratios")

rr_uf <- rbindlist(lapply(UF, function(u) rbind(
  rr_row("State", u, by_year(cn[uf == u]), pop_uf[, u]),
  rr_row("State, excluding PDL", u, by_year(cn[uf == u & (is.na(pdl) | pdl != 1)]), pop_uf[, u]))))
save_table(rr_uf, "rate_ratios_states")

# Rate ratios by race/skin colour versus White, whole period, with interpolated and with
# fixed 2022 denominators (Supplementary Table S5)
race_tab <- rbindlist(lapply(setdiff(RACES, "Asian"), function(rc) {
  n <- by_year(cn[race == rc])
  data.table(race = rc, n = sum(n), rate = sum(n) / sum(pop_race[, rc]) * 1e5,
             rate_fixed2022 = sum(n) / (pop_race_fixed[rc] * length(YEARS)) * 1e5)
}))
nW <- race_tab[race == "White", n]
race_tab[, se := sqrt(1 / n + 1 / nW)]
race_tab[, `:=`(RR_vs_white = rate / rate[race == "White"], RR_vs_white_fixed2022 = rate_fixed2022 / rate_fixed2022[race == "White"])]
race_tab[, `:=`(RR_lo = RR_vs_white * exp(-1.96 * se), RR_hi = RR_vs_white * exp(1.96 * se),
                RR_fixed_lo = RR_vs_white_fixed2022 * exp(-1.96 * se), RR_fixed_hi = RR_vs_white_fixed2022 * exp(1.96 * se))]
race_tab <- merge(race_tab, rr[group == "Race/skin colour", .(race = category, rate_Pre, rate_Post, RR_Post, RR_Post_lo, RR_Post_hi)], by = "race")
# Absolute inequality: rate difference versus White (per 100,000 per year) before and after the pandemic
race_tab[, `:=`(RD_vs_white_Pre = rate_Pre - rate_Pre[race == "White"], RD_vs_white_Post = rate_Post - rate_Post[race == "White"])]
race_tab[, pct_missing_race_overall := mean(is.na(cn$race)) * 100]
race_tab <- merge(race_tab, rbindlist(lapply(setdiff(RACES, "Asian"), function(rc) {
  r <- rr_row("Race/skin colour", rc, by_year(cn[race == rc]), pop_race_extrap[, rc])
  data.table(race = rc, RR_Post_extrap = r$RR_Post, RR_Post_extrap_lo = r$RR_Post_lo, RR_Post_extrap_hi = r$RR_Post_hi)
})), by = "race")
save_table(race_tab[, !"se"], "table_s5_race")

# Ascertainment check by age: rate ratios (2022-2025 vs 2015-2019) for all new pulmonary cases,
# with or without laboratory confirmation, and the share of new pulmonary cases that were confirmed.
# Rapid molecular testing raised confirmation more in children than in adults, so the increase in
# laboratory-confirmed cases among children overstates the increase in diagnosed disease.
new_pulm <- all_pulm[TRATAMENTO == "1"]
new_pulm[, period := fifelse(year <= 2019, "Pre", fifelse(year <= 2021, "Pandemic", "Post"))]
age_check <- rbindlist(lapply(c("All ages", AGES), function(a) {
  x <- if (a == "All ages") new_pulm else new_pulm[age_group == a]
  pa <- if (a == "All ages") pop_br else pop_age_br[, a]
  lab_rr <- rr_row("lab", a, by_year(x[lab == TRUE]), pa)
  all_rr <- rr_row("all", a, by_year(x), pa)
  data.table(age_group = a,
             pct_confirmed_Pre = x[period == "Pre", mean(lab) * 100], pct_confirmed_Post = x[period == "Post", mean(lab) * 100],
             RR_lab_confirmed = lab_rr$RR_Post, RR_lab_lo = lab_rr$RR_Post_lo, RR_lab_hi = lab_rr$RR_Post_hi,
             RR_all_new_pulmonary = all_rr$RR_Post, RR_all_lo = all_rr$RR_Post_lo, RR_all_hi = all_rr$RR_Post_hi)
}))
save_table(age_check, "rate_ratios_by_age_case_definition")

# ---- 8. Characteristics, laboratory indicators and HIV by period --------------

describe <- function(var, levels) {
  tab <- table(factor(fifelse(is.na(cn[[var]]), "Missing", as.character(cn[[var]])),
                      levels = c(levels, "Missing")), factor(cn$period, levels = names(PERIOD)))
  known <- tab[rownames(tab) != "Missing", , drop = FALSE]; known <- known[rowSums(known) > 0, , drop = FALSE]
  ct <- suppressWarnings(chisq.test(known))
  V <- sqrt(ct$statistic / (sum(known) * (min(dim(known)) - 1)))
  out <- as.data.table(as.data.frame.matrix(tab), keep.rownames = "level")
  for (p in names(PERIOD)) out[, (paste0(p, "_pct")) := get(p) / sum(get(p)) * 100]
  cbind(variable = var, out, p_value = ct$p.value, cramers_v = as.numeric(V))
}
cn[, `:=`(pdl_yn = fcase(pdl == 1, "Yes", pdl == 0, "No"), homeless_yn = fcase(homeless == 1, "Yes", homeless == 0, "No"))]
table1 <- rbindlist(list(describe("sex", c("Male", "Female")), describe("age_group", AGES),
                         describe("race", c("White", "Black", "Brown (pardo)", "Asian", "Indigenous")),
                         describe("region", REGION), describe("pdl_yn", c("Yes", "No")),
                         describe("homeless_yn", c("Yes", "No"))))
save_table(table1, "table1_characteristics")
comorbidities <- cn[, lapply(.SD, function(v) mean(v, na.rm = TRUE) * 100), by = period,
                    .SDcols = c("alcohol", "drugs", "tobacco", "diabetes")]
save_table(comorbidities, "comorbidities_by_period")

dst_done <- c("1", "2", "3", "4", "5")
lab_hiv <- cn[, {
  rbind(cbind(indicator = "Rapid molecular test done", wilson(sum(TEST_MOLEC %chin% c("1", "2", "3", "4")), .N)),
        cbind(indicator = "Culture done", wilson(sum(CULTURA_ES %chin% c("1", "2", "3")), .N)),
        cbind(indicator = "Drug susceptibility test among culture-positive", wilson(sum(CULTURA_ES == "1" & TEST_SENSI %chin% dst_done, na.rm = TRUE), sum(CULTURA_ES %chin% "1"))),
        cbind(indicator = "Rifampicin resistance among tested (DST)", wilson(sum(TEST_SENSI %chin% c("2", "3")), sum(TEST_SENSI %chin% dst_done))),
        cbind(indicator = "Multidrug resistance among tested (DST)", wilson(sum(TEST_SENSI %chin% "3"), sum(TEST_SENSI %chin% dst_done))),
        cbind(indicator = "Rifampicin resistance on rapid molecular test", wilson(sum(TEST_MOLEC %chin% "2"), sum(TEST_MOLEC %chin% c("1", "2")))),
        cbind(indicator = "HIV test result available", wilson(sum(hiv %chin% c("Positive", "Negative")), .N)),
        cbind(indicator = "HIV positive among tested", wilson(sum(hiv %chin% "Positive"), sum(hiv %chin% c("Positive", "Negative")))))
}, by = period]
save_table(lab_hiv, "laboratory_hiv_by_period")

# ---- 9. Treatment outcomes and loss to follow-up ----------------------------

# Outcome codes (SITUA_ENCE): 1 cure, 2 loss to follow-up, 3 death from TB, 4 death from other
# causes, 5 transfer, 6 change of diagnosis, 7 drug-resistant TB, 8 change of regimen, 9 failure,
# 10 primary loss to follow-up. Any code 1-10 counts as a closed record.
cn[, `:=`(closed = !is.na(outcome) & outcome %between% c(1, 10), ltfu = outcome %in% c(2, 10))]
outcomes_year <- cn[, .(n = .N, pct_not_closed = mean(!closed) * 100,
                        cure = sum(outcome == 1, na.rm = TRUE) / sum(closed) * 100,
                        ltfu = sum(ltfu) / sum(closed) * 100,
                        death_tb = sum(outcome == 3, na.rm = TRUE) / sum(closed) * 100,
                        death_other = sum(outcome == 4, na.rm = TRUE) / sum(closed) * 100,
                        transfer = sum(outcome == 5, na.rm = TRUE) / sum(closed) * 100), keyby = year]
save_table(outcomes_year, "outcomes_by_year")
outcome_apc <- rbind(cbind(indicator = "Cure", prais_winsten_apc(2015:2024, outcomes_year[year <= 2024, cure])),
                     cbind(indicator = "Loss to follow-up", prais_winsten_apc(2015:2024, outcomes_year[year <= 2024, ltfu])))
save_table(outcome_apc, "outcomes_apc")
# Sensitivity: loss to follow-up with all cases as denominator (records still open counted as not lost)
ltfu_all_cases <- cn[year <= 2024, .(ltfu_pct_closed = sum(ltfu) / sum(closed) * 100, ltfu_pct_all = mean(ltfu) * 100), keyby = year]
save_table(ltfu_all_cases, "ltfu_denominator_sensitivity")
save_table(cbind(indicator = "Loss to follow-up, all cases as denominator",
                 prais_winsten_apc(2015:2024, ltfu_all_cases$ltfu_pct_all)), "outcomes_apc_all_cases")

# Closed cohorts diagnosed up to 2024 (2025 records are mostly still open)
coh <- cn[year <= 2024 & closed == TRUE]
outcomes_period <- coh[, rbind(cbind(outcome = "Cure", wilson(sum(outcome == 1), .N)),
                               cbind(outcome = "Loss to follow-up", wilson(sum(ltfu), .N)),
                               cbind(outcome = "Death from TB", wilson(sum(outcome == 3), .N))), by = period]
subgroups <- list(PDL = coh$pdl %in% 1, Homeless = coh$homeless %in% 1, "HIV positive" = coh$hiv %chin% "Positive",
                  "Illicit drugs" = coh$drugs %in% 1,
                  Others = !(coh$pdl %in% 1) & !(coh$homeless %in% 1) & !(coh$hiv %chin% "Positive"))
outcomes_subgroups <- rbindlist(lapply(subgroups, function(k) rbind(
  cbind(outcome = "Cure", n = sum(k), wilson(sum(coh$outcome[k] == 1), sum(k))),
  cbind(outcome = "Loss to follow-up", n = sum(k), wilson(sum(coh$ltfu[k]), sum(k))))), idcol = "subgroup")
save_table(outcomes_period, "outcomes_by_period"); save_table(outcomes_subgroups, "outcomes_by_subgroup")
save_table(coh[, .(n = .N, cure = mean(outcome == 1) * 100, ltfu = mean(ltfu) * 100), by = .(uf, period)], "outcomes_by_state_period")

# Modified Poisson regression (log link, robust variance) for loss to follow-up, transfers excluded,
# state fixed effects and standard errors clustered by state. Missing covariates coded as "No"
# (binary) or as a separate category, as in the manuscript.
rg <- coh[outcome != 5 & !is.na(sex)]
rg[, `:=`(y = as.integer(ltfu),
          sex_f = relevel(factor(sex), "Female"),
          age_f = relevel(factor(fifelse(is.na(age_group), "25-34", fifelse(age_group %chin% c("0-14", "15-24"), "0-24", age_group))), "35-44"),
          race_f = relevel(factor(fcase(is.na(race), "Missing", race %chin% c("Asian", "Indigenous"), "Asian/Indigenous", rep(TRUE, .N), race)), "White"),
          hiv_f = relevel(factor(fifelse(is.na(hiv) | hiv %chin% c("Not done", "Pending"), "Not done/missing", hiv)), "Negative"),
          dot_f = relevel(factor(fcase(TRAT_SUPER == "1", "Yes", TRAT_SUPER == "2", "No", default = "Missing")), "Yes"),
          period_f = factor(period, levels = names(PERIOD)), uf_f = factor(uf))]
for (v in c("pdl", "homeless", "alcohol", "drugs", "tobacco", "diabetes")) rg[, (paste0(v, "_b")) := fifelse(is.na(get(v)), 0L, get(v))]
ltfu_fit <- glm(y ~ sex_f + age_f + race_f + pdl_b + homeless_b + hiv_f + alcohol_b + drugs_b + tobacco_b + diabetes_b +
                  dot_f + period_f + uf_f, family = poisson(), data = rg)
V_cl <- vcovCL(ltfu_fit, cluster = ~uf_f, type = "HC0")
ltfu_ct <- coeftest(ltfu_fit, vcov. = V_cl)
ltfu_tab <- data.table(term = rownames(ltfu_ct), PR = exp(ltfu_ct[, 1]), lo = exp(ltfu_ct[, 1] - 1.96 * ltfu_ct[, 2]),
                       hi = exp(ltfu_ct[, 1] + 1.96 * ltfu_ct[, 2]), p = ltfu_ct[, 4])[!grepl("^uf_f|Intercept", term)]
ltfu_tab[, `:=`(n = nrow(rg), prevalence_pct = mean(rg$y) * 100)]
save_table(ltfu_tab, "table_s4_ltfu_model")

# ---- 10. Municipal and spatial analysis --------------------------------------

# Municipal counts (2015-2025) against the 2022 census population, multiplied by the number of years.
mun <- fread("data/municipal_population_census2022.csv", colClasses = list(character = "cod7"))
mun[, mun6 := substr(cod7, 1, 6)]
count_mun <- function(x) { v <- x[, .N, by = mun6]; out <- v$N[match(mun$mun6, v$mun6)]; out[is.na(out)] <- 0; out }
nonpdl <- cn[is.na(pdl) | pdl != 1]
mun[, `:=`(cases = count_mun(cn), cases_pdl = count_mun(cn[pdl == 1]), cases_np = count_mun(nonpdl),
           cases_pre_np = count_mun(nonpdl[year <= 2019]), cases_post_np = count_mun(nonpdl[year >= 2022]))]
mun[, `:=`(rate_np = cases_np / (pop2022 * 11) * 1e5, pct_pdl = fifelse(cases > 0, cases_pdl / cases * 100, NA_real_))]
unmatched_cases <- cn[!mun6 %chin% mun$mun6, .N]

# Geometry: municipal boundaries from the geodata-br project (IBGE 2010 meshes, 5,564 polygons
# that match the 2022 municipal codes used here except for five municipalities created after 2010).
GEO_URL  <- "https://raw.githubusercontent.com/tbrugz/geodata-br/master/geojson/geojs-100-mun.json"
GEO_FILE <- "data-raw/geojs-100-mun.json"
if (!file.exists(GEO_FILE)) download.file(GEO_URL, GEO_FILE, mode = "wb")
geo <- st_read(GEO_FILE, quiet = TRUE)
geo <- geo[, c("id", "name")]; names(geo)[1] <- "cod7"
geo <- st_make_valid(geo)
g <- merge(geo, mun, by = "cod7")                       # keeps municipalities with geometry and population
g <- g[order(g$cod7), ]

# Queen contiguity; municipalities without neighbours (islands) get their nearest neighbour.
sf_use_s2(FALSE)
nb <- poly2nb(g, queen = TRUE, snap = 1e-6)
islands <- which(card(nb) == 0)
if (length(islands)) {
  knn1 <- knn2nb(knearneigh(st_coordinates(suppressWarnings(st_point_on_surface(st_geometry(g)))), k = 1))
  for (i in islands) { nb[[i]] <- knn1[[i]]; j <- knn1[[i]]; nb[[j]] <- sort(unique(c(setdiff(nb[[j]], 0L), i))) }
}
lw <- nb2listw(nb, style = "W")

uf_of <- function(cod) unname(UF[substr(cod, 1, 2)])
person_years <- g$pop2022 * length(YEARS)
# Empirical Bayes smoothing (global and spatial, Marshall 1991)
g$eb_np  <- EBest(g$cases_np, person_years)$estmm * 1e5
g$seb_np <- EBlocal(g$cases_np, person_years, nb)$est * 1e5
g$eb_all <- EBest(g$cases, person_years)$estmm * 1e5

# Moran's I for rates (Assuncao-Reis empirical Bayes standardisation), 999 permutations
moran_np  <- EBImoran.mc(g$cases_np, person_years, lw, nsim = N_PERM)
moran_all <- EBImoran.mc(g$cases, person_years, lw, nsim = N_PERM)

# Local Moran's I on Assuncao-Reis standardised rates, conditional permutation.
ar_z <- function(cases, pop) {
  b <- sum(cases) / sum(pop); r <- cases / pop
  s2 <- sum(pop * (r - b)^2) / sum(pop); a <- s2 - b / mean(pop)
  v <- a + b / pop; v[v < 0] <- (b / pop)[v < 0]
  (r - b) / sqrt(v)
}
lisa <- function(cases, nsim) {
  z <- ar_z(cases, person_years)
  lm <- localmoran_perm(z, lw, nsim = nsim, iseed = SEED, alternative = "two.sided")
  q <- attr(lm, "quadr")$mean
  list(p = lm[, "Pr(folded) Sim"], quadrant = as.character(q))
}
l999 <- lisa(g$cases_np, N_PERM)
g$lisa_np <- ifelse(l999$p < 0.05, l999$quadrant, "Not significant")
l9999 <- lisa(g$cases_np, 9999)
hh_raw <- l9999$quadrant == "High-High" & l9999$p < 0.05
hh_fdr <- l9999$quadrant == "High-High" & p.adjust(l9999$p, "BH") < 0.05
g$lisa_np_fdr <- ifelse(hh_fdr, "High-High (FDR)", "Not significant")
spatial_summary <- data.table(
  measure = c("Municipalities analysed", "Cases without municipal match", "Global Moran's I (EB), non-PDL", "Permutation p",
              "Global Moran's I (EB), all cases", "Permutation p",
              "High-high municipalities, p<0.05 (999 perm.)",
              "High-high municipalities, p<0.05 (9,999 perm.)", "  share of population (%)", "  share of non-PDL cases (%)",
              "High-high municipalities, FDR<0.05", "  share of population (%)", "  share of non-PDL cases (%)"),
  value = c(nrow(g), unmatched_cases, moran_np$statistic, moran_np$p.value, moran_all$statistic, moran_all$p.value,
            sum(g$lisa_np == "High-High"),
            sum(hh_raw), sum(g$pop2022[hh_raw]) / sum(g$pop2022) * 100, sum(g$cases_np[hh_raw]) / sum(g$cases_np) * 100,
            sum(hh_fdr), sum(g$pop2022[hh_fdr]) / sum(g$pop2022) * 100, sum(g$cases_np[hh_fdr]) / sum(g$cases_np) * 100))
save_table(spatial_summary, "spatial_summary")
lisa_counts <- as.data.table(st_drop_geometry(g))[, .(municipalities = .N, population = sum(pop2022), cases_np = sum(cases_np)),
                                                   keyby = .(lisa_np, uf = uf_of(cod7))]
save_table(lisa_counts, "lisa_clusters_by_state")
save_table(as.data.table(st_drop_geometry(g))[, .(cod7, name, rgi, pop2022, cases, cases_np, cases_pdl, cases_pre_np, cases_post_np,
                                                  pct_pdl, rate_np, eb_all, eb_np, seb_np, lisa_np, lisa_np_fdr)], "municipal_results")

# Municipalities with >= 100,000 inhabitants and the highest non-PDL rates (Table S6),
# with the post- versus pre-pandemic rate ratio
top_mun <- as.data.table(st_drop_geometry(g))[pop2022 >= 1e5][order(-rate_np)][1:25]
top_mun[, `:=`(uf = uf_of(cod7), RR = (cases_post_np / 4) / (cases_pre_np / 5))]
top_mun[, `:=`(RR_lo = RR * exp(-1.96 * sqrt(1 / cases_post_np + 1 / cases_pre_np)),
               RR_hi = RR * exp(1.96 * sqrt(1 / cases_post_np + 1 / cases_pre_np)))]
save_table(top_mun[, .(cod7, name, uf, pop2022, cases_np, rate_np, eb_np, pct_pdl, RR, RR_lo, RR_hi, lisa_np)], "table_s6_top_municipalities")

# Funnel plot control limits. Poisson limits ignore the extra variation between municipalities, so
# far more units fall outside them than chance would allow. Limits adjusted for overdispersion use the
# additive random-effects method of Spiegelhalter (2005): z-scores winsorised at the 10th and 90th
# percentiles give the overdispersion factor, from which the between-municipality variance tau^2 follows.
funnel_limits <- function(cases, py) {
  r0 <- sum(cases) / sum(py); r <- cases / py; s <- sqrt(r0 / py)
  z <- (r - r0) / s; q <- quantile(z, c(0.1, 0.9)); zw <- pmin(pmax(z, q[1]), q[2])
  phi <- mean(zw^2); w <- 1 / s^2; N <- length(z)
  tau2 <- max(0, (N * phi - (N - 1)) / (sum(w) - sum(w^2) / sum(w)))
  list(r0 = r0, tau2 = tau2, phi = phi, s = s, r = r)
}
fl <- funnel_limits(g$cases_np, person_years)
big <- g$pop2022 >= 1e5
funnel_summary <- data.table(
  measure = c("Municipalities >= 100,000 inhabitants", "Above Poisson 99.8% limit", "Above overdispersion-adjusted 99.8% limit",
              "Below overdispersion-adjusted 99.8% limit", "Overdispersion factor (winsorised)", "tau (per 100,000)"),
  value = c(sum(big), sum(big & fl$r > fl$r0 + 3.09 * fl$s), sum(big & fl$r > fl$r0 + 3.09 * sqrt(fl$s^2 + fl$tau2)),
            sum(big & fl$r < fl$r0 - 3.09 * sqrt(fl$s^2 + fl$tau2)), fl$phi, sqrt(fl$tau2) * 1e5))
save_table(funnel_summary, "funnel_summary")
g$above_adjusted_998 <- fl$r > fl$r0 + 3.09 * sqrt(fl$s^2 + fl$tau2)

# Immediate geographic regions (IBGE 2017): non-PDL rates and post/pre rate ratios
rgi_names <- fread("data/immediate_regions.csv")
rgi <- mun[, .(pop2022 = sum(pop2022), cases_np = sum(cases_np), cases_pre_np = sum(cases_pre_np),
               cases_post_np = sum(cases_post_np), cases_pdl = sum(cases_pdl)), by = rgi]
rgi[, rate_np := cases_np / (pop2022 * 11) * 1e5]
rgi[, c("RR_np", "RR_np_lo", "RR_np_hi") := rate_ratio(pmax(cases_post_np, 0.5), pop2022 * 4, pmax(cases_pre_np, 0.5), pop2022 * 5)]
rgi <- merge(rgi, rgi_names, by = "rgi")
save_table(rgi[order(-rate_np)], "immediate_regions")

# ---- 11. Tuberculosis mortality (SIM) -----------------------------------------

# Deaths with tuberculosis (ICD-10 A15-A19) as the underlying cause, by month of death and state
# of residence, from the DATASUS TabNet query described in data-raw/README.md.
cf_deaths <- counterfactual(mf$deaths)
deaths_sens <- rbindlist(list(
  "Main" = cf_deaths$summary,
  "No trend, baseline Mar 2017-Feb 2020" = counterfactual(mf$deaths, trend = FALSE, start = 26)$summary,
  "No trend, baseline Mar 2018-Feb 2020" = counterfactual(mf$deaths, trend = FALSE, start = 38)$summary,
  "Post-pandemic window 2022-2024" = counterfactual(mf$deaths, end = 120)$summary), idcol = "specification")
save_table(deaths_sens, "deaths_counterfactual")
new_all_forms <- by_year(d[year %between% c(2015, 2025) & TRATAMENTO == "1" & !is.na(ID_MN_RESI)])
mortality <- data.table(year = YEARS, deaths = deaths_monthly[, sum(BR), keyby = year]$V1)
mortality[, `:=`(mortality_rate = deaths / pop_br * 1e5, new_cases_all_forms = new_all_forms,
                 deaths_per_100_new_cases = deaths / new_all_forms * 100)]
save_table(mortality, "mortality_annual")
save_table(fread("data/sim_tb_deaths_by_icd_year.csv", header = TRUE), "deaths_by_icd_year")

# Population ageing. Deaths by age group and year (TabNet, same query with age as columns) give
# (i) age-standardised mortality rates (direct method, 2022 population; deaths of unknown age
# distributed in proportion to known ages) and (ii) a counterfactual whose offset is the number of
# deaths expected from the age structure of each year at the 2015-2019 age-specific rates, so that
# expected deaths grow with ageing and not only with population size.
deaths_age <- fread("data/sim_tb_deaths_by_age_year.csv", header = TRUE)
da <- as.matrix(deaths_age[match(YEARS, year), AGES, with = FALSE]); rownames(da) <- YEARS
da <- da * (deaths_age$total / rowSums(da))                                    # redistribute unknown age
mortality[, asmr := as.vector((da / pop_age_br * 1e5) %*% std_weights)]
ref_rates <- colSums(da[as.character(2015:2019), ]) / colSums(pop_age_br[as.character(2015:2019), ])
expected_by_age <- as.vector(pop_age_br %*% ref_rates)
pop_age_adjusted <- expected_by_age / (sum(da[as.character(2015:2019), ]) / sum(pop_br[as.character(2015:2019)]))
mf$deaths_age_adj <- monthly_frame(deaths_monthly$BR, setNames(pop_age_adjusted, YEARS))
cf_deaths_age <- counterfactual(mf$deaths_age_adj)
deaths_sens <- rbind(deaths_sens, cbind(specification = "Offset adjusted for population ageing", cf_deaths_age$summary))
save_table(deaths_sens, "deaths_counterfactual")
save_table(mortality, "mortality_annual")
mortality_apc <- rbind(cbind(series = "Crude mortality", window = "2015-2019", prais_winsten_apc(2015:2019, mortality[year <= 2019, mortality_rate])),
                       cbind(series = "Age-standardised mortality", window = "2015-2019", prais_winsten_apc(2015:2019, mortality[year <= 2019, asmr])))
save_table(mortality_apc, "mortality_apc")

# ---- 12. Figures ---------------------------------------------------------------

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; INK <- "#52514e"; GREY <- "#8a8984"
theme_set(theme_minimal(base_size = 9) +
            theme(panel.grid.minor = element_blank(), plot.title.position = "plot",
                  plot.title = element_text(face = "bold", size = 10), strip.text = element_text(face = "bold", hjust = 0)))
save_fig <- function(p, name, w = 7.4, h = 5) ggsave(file.path("output/figures", paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")

# Figure 1: observed and expected monthly counts with 95% prediction intervals
series_label <- c(all = "All new laboratory-confirmed pulmonary cases", non_pdl = "Excluding persons deprived of liberty",
                  pdl = "Persons deprived of liberty (PDL)", deaths = "Deaths with TB as underlying cause (SIM)")
res_cf <- c(cf, list(deaths = cf_deaths))
fig1_data <- rbindlist(lapply(names(series_label), function(k)
  data.table(series = series_label[k], date = as.Date(months), observed = mf[[k]]$y,
             expected = res_cf[[k]]$expected, lo = res_cf[[k]]$pi_lo, hi = res_cf[[k]]$pi_hi)))
fig1_data[, series := factor(series, levels = series_label)]
fig1 <- ggplot(fig1_data, aes(date)) +
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2022-01-01"), ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = GREY, alpha = 0.25) +
  geom_line(aes(y = expected, linetype = "Expected (pre-pandemic model)", colour = "Expected (pre-pandemic model)"), linewidth = 0.5) +
  geom_line(aes(y = observed, linetype = "Observed", colour = "Observed"), linewidth = 0.5) +
  scale_linetype_manual(NULL, values = c("Observed" = "solid", "Expected (pre-pandemic model)" = "22")) +
  scale_colour_manual(NULL, values = c("Observed" = BLUE, "Expected (pre-pandemic model)" = INK)) +
  facet_wrap(~series, ncol = 2, scales = "free_y") +
  labs(x = NULL, y = "Monthly count", title = "Observed and expected tuberculosis notifications and deaths, Brazil, 2015-2025",
       caption = "Grey band: pandemic period (March 2020-December 2021). Shaded ribbon: 95% prediction interval.") +
  theme(legend.position = "top")
save_fig(fig1, "figure1_observed_expected", h = 5.4)

# Figure 2: annual crude and age-standardised rates, and rates by region
fig2_data <- rbind(annual[, .(year, series = "Brazil, crude", rate)], annual[, .(year, series = "Brazil, age-standardised", rate = asr)],
                   annual[, .(year, series = "Brazil, excluding PDL", rate = rate_non_pdl)],
                   rbindlist(lapply(REGION, function(r) data.table(year = YEARS, series = r, rate = region_rate[, r]))))
fig2_data[, panel := fifelse(grepl("^Brazil", series), "A. Brazil", "B. Regions")]
fig2 <- ggplot(fig2_data, aes(year, rate, colour = series)) +
  annotate("rect", xmin = 2019.5, xmax = 2021.5, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_line(linewidth = 0.6) + geom_point(size = 1) + facet_wrap(~panel, scales = "free_y") +
  scale_x_continuous(breaks = seq(2015, 2025, 2)) +
  labs(x = NULL, y = "Rate per 100,000", colour = NULL,
       title = "New laboratory-confirmed pulmonary tuberculosis, rates per 100,000, Brazil, 2015-2025") +
  theme(legend.position = "bottom")
save_fig(fig2, "figure2_annual_rates", h = 4.2)

# Figure 3: post- versus pre-pandemic rate ratios
fig3_data <- rr[group != "Total" | category == "Brazil"][, label := factor(category, levels = rev(category))]
fig3 <- ggplot(fig3_data, aes(RR_Post, label, colour = RR_Post > 1)) +
  geom_vline(xintercept = 1, colour = GREY) +
  geom_pointrange(aes(xmin = RR_Post_lo, xmax = RR_Post_hi), size = 0.25) +
  scale_colour_manual(values = c(`TRUE` = ORANGE, `FALSE` = BLUE), guide = "none") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Rate ratio, 2022-2025 vs 2015-2019 (95% CI)", y = NULL,
       title = "Change in notification rates after the pandemic, by group") +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0))
save_fig(fig3, "figure3_rate_ratios", w = 6, h = 5.5)

# Figure 4: smoothed municipal rates and local clusters (non-PDL cases)
gf <- suppressWarnings(st_simplify(g, dTolerance = 0.01))   # for drawing only
brks <- unique(quantile(gf$seb_np, probs = seq(0, 1, 0.2), na.rm = TRUE))
gf$seb_class <- cut(gf$seb_np, brks, include.lowest = TRUE, dig.lab = 4)
fig4a <- ggplot(gf) + geom_sf(aes(fill = seb_class), colour = NA) +
  scale_fill_brewer(palette = "OrRd", name = "Spatial EB rate\nper 100,000") +
  labs(title = "A. Smoothed rate, excluding PDL") + theme_void(base_size = 9)
gf$lisa_plot <- factor(fifelse(gf$lisa_np == "Not significant", "Not significant", gf$lisa_np),
                       levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant"))
fig4b <- ggplot(gf) + geom_sf(aes(fill = lisa_plot), colour = NA) +
  scale_fill_manual(values = c("High-High" = "#c0392b", "Low-Low" = BLUE, "High-Low" = "#f1a07c",
                               "Low-High" = "#9cc3ef", "Not significant" = "grey90"), name = "Local Moran") +
  labs(title = "B. Local clusters (p < 0.05, 999 permutations)") + theme_void(base_size = 9)
save_fig(fig4a, "figure4a_eb_rate_map", w = 5, h = 5); save_fig(fig4b, "figure4b_lisa_map", w = 5, h = 5)

# Supplementary figure: funnel plot of municipal non-PDL rates
br_rate <- sum(g$cases_np) / sum(person_years)
pp <- 10^seq(log10(min(g$pop2022)), log10(max(g$pop2022)), length.out = 300)
limits <- rbindlist(list(
  "Poisson 99.8%" = data.table(pop = pp, sd = sqrt(br_rate / (pp * 11))),
  "Overdispersion-adjusted 99.8%" = data.table(pop = pp, sd = sqrt(br_rate / (pp * 11) + fl$tau2))), idcol = "limit")
limits[, `:=`(lo = pmax(br_rate - 3.09 * sd, 0) * 1e5, hi = (br_rate + 3.09 * sd) * 1e5)]
figs <- ggplot() +
  geom_point(data = st_drop_geometry(g), aes(pop2022, rate_np, colour = pop2022 >= 1e5), size = 0.6, alpha = 0.5) +
  geom_line(data = limits, aes(pop, hi, linetype = limit), colour = INK) +
  geom_line(data = limits, aes(pop, lo, linetype = limit), colour = INK) +
  geom_hline(yintercept = br_rate * 1e5, colour = ORANGE) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = c(`TRUE` = BLUE, `FALSE` = GREY), labels = c("< 100,000", ">= 100,000"), name = "Population") +
  labs(x = "Population, 2022 census (log scale)", y = "Rate per 100,000 person-years, excluding PDL", linetype = "Control limits",
       title = "Municipal notification rates, 2015-2025, funnel plot")
save_fig(figs, "figure_s4_funnel", w = 6.5, h = 4.5)

# ---- 13. Save results and session information ----------------------------------

results <- list(flow = flow, validation = validation, annual = annual, apc = apc, its = its_tab, its_trends = its_trends,
                its_levels = its_levels, counterfactual = cf_table, counterfactual_geo = cf_geo, sensitivity = sens_table,
                autocorrelation = autocorr, rate_ratios = rr, rate_ratios_states = rr_uf, race = race_tab, table1 = table1,
                laboratory_hiv = lab_hiv, outcomes_year = outcomes_year, outcomes_period = outcomes_period,
                outcomes_subgroups = outcomes_subgroups, ltfu_model = ltfu_tab, spatial = spatial_summary,
                top_municipalities = top_mun, immediate_regions = rgi, mortality = mortality, deaths_counterfactual = deaths_sens)
saveRDS(results, "output/results.rds")
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")

# Key numbers quoted in the manuscript summary
cat("\nNew laboratory-confirmed pulmonary cases:", nrow(cn), "\n")
print(cf_table[window %in% c("Pandemic", "Post", "Cumulative")][, .(series, window, observed, expected = round(expected), OE = round(OE, 3), OE_lo = round(OE_lo, 3), OE_hi = round(OE_hi, 3))])
print(cf_deaths$summary[window %in% c("Pandemic", "Post", "Cumulative"), .(window, OE = round(OE, 3), OE_lo = round(OE_lo, 3), OE_hi = round(OE_hi, 3))])
print(race_tab[, .(race, RR_vs_white = round(RR_vs_white, 2), RR_Post = round(RR_Post, 2))])
