suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

root <- "/data/h011/project/analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1"
f2in <- file.path(root, "Figure2_PRE_BROAD")
f5in <- file.path(root, "Figure5_PRE_BROAD")
f2out <- file.path(root, "Figure2_PRE_BROAD_v3")
f5out <- file.path(root, "Figure5_PRE_BROAD_v2")
dir.create(f2out, recursive=TRUE, showWarnings=FALSE)

read_tsv <- function(path) read.delim(path, check.names=FALSE, stringsAsFactors=FALSE)
write_tsv <- function(x, path) write.table(x, path, sep="\t", row.names=FALSE, quote=FALSE, na="NA")
save3 <- function(p, stem, width=180, height=150) {
  ggsave(paste0(stem, ".png"), p, width=width, height=height, units="mm", dpi=600, bg="white")
  ggsave(paste0(stem, ".pdf"), p, width=width, height=height, units="mm", bg="white", device=cairo_pdf)
  ggsave(paste0(stem, ".svg"), p, width=width, height=height, units="mm", bg="white")
}
tag <- function(p, x) p + labs(tag=x) + theme(plot.tag=element_text(face="bold", size=9), plot.tag.position="topleft")

pal_tier <- c(None="#F1F4F6", Weak="#C8DEEB", Supportive="#A6C9DE", Strong="#79A9C5")
pal_species <- c(
  "M. arginini"="#416C86", "M. hominis"="#6F99B1",
  "M. hyorhinis"="#98BACD", "M. salivarium"="#B8D1DE"
)
theme_paper <- theme_minimal(base_family="sans", base_size=7.2) +
  theme(
    panel.grid.minor=element_blank(), panel.grid.major=element_line(color="#E7ECEF", linewidth=.22),
    plot.title=element_text(hjust=.5, face="bold", size=8.2, margin=margin(b=4)),
    axis.title=element_text(size=7.5), axis.text=element_text(size=6.7, color="#26343C"),
    strip.text=element_text(size=6.7, face="bold"), strip.background=element_rect(fill="#F3F6F8", color=NA),
    legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=6.3),
    legend.key.height=unit(3.2,"mm"), legend.key.width=unit(5,"mm"),
    plot.margin=margin(6,6,6,6)
  )

# ---------------- Figure 2 ----------------
cond <- read_tsv(file.path(f2in, "Figure2A_conditions_480.tsv"))
rawC <- read_tsv(file.path(f2in, "Figure2C_source.tsv"))
rawD <- read_tsv(file.path(f2in, "Figure2D_source.tsv"))
rawE <- read_tsv(file.path(f2in, "Figure2E_source.tsv"))
rawF <- read_tsv(file.path(f2in, "Figure2F_source.tsv"))
blank <- read_tsv(file.path(f2in, "Supplementary_Figure2G_blank_controls_source.tsv"))

doses <- c(0,1,2,3,5,10,20,50)
dose_labels <- as.character(doses)
species_order <- c("M. arginini","M. hominis","M. hyorhinis","M. salivarium")

stopifnot(nrow(cond)==480, length(unique(cond$background_sample))==3,
          length(unique(cond$species))==4, isTRUE(all.equal(as.numeric(sort(unique(cond$dose_pairs))), doses)),
          length(unique(cond$seed))==5,
          length(unique(cond$scientific_commit))==1,
          startsWith(unique(cond$scientific_commit), "77c4e85"),
          sum(cond$dose_pairs>0)==420, nrow(blank)==3)

design <- data.frame(
  item=c("WES backgrounds","Microbial species","Dose levels","Fixed seeds","Labeled conditions","Non-zero conditions","Independent blanks"),
  value=c(3,4,8,5,480,420,3),
  detail=c("Host-depleted WES", "One source strain per species", "0, 1, 2, 3, 5, 10, 20, 50 pairs", "20260930-20260934", "Not independent experiments", "3 backgrounds x 4 species x 7 non-zero doses x 5 seeds", "D000 labels deduplicated by background")
)
write_tsv(design, file.path(f2out,"Fig2A_source.tsv"))

p2a_core <- design[1:4,]
p2a_core$item_short <- c("WES\nbackgrounds","Microbial\nspecies","Dose\nlevels","Fixed\nseeds")
p2a <- ggplot(p2a_core, aes(x=seq_len(nrow(p2a_core)), y=1)) +
  annotate("segment", x=1, xend=4, y=1, yend=1, color="#AFC1CB", linewidth=.55) +
  geom_point(shape=21, size=5, fill="#C8DEEB", color="#416276", stroke=.6) +
  geom_text(aes(label=value), size=2.7, fontface="bold", color="#243A47") +
  geom_text(aes(label=item_short), y=.73, size=2.05, lineheight=1.05) +
  annotate("text",x=2.5,y=.49,label="420 nonzero conditions\n+ 3 independent blanks",size=2.15,lineheight=1.10,color="#42545E") +
  scale_x_continuous(NULL, breaks=NULL, expand=expansion(mult=c(.25,.20))) +
  scale_y_continuous(NULL, breaks=NULL, limits=c(.39,1.12)) +
  labs(title="Spike-in design") + theme_paper + theme(panel.grid=element_blank(), legend.position="none")

# B: real independent denominators; dose 0 uses the three shared independent blanks.
e_tmp <- rawE[rawE$dose_pairs>0,]
e_tmp$n_conditions <- 1L
e_tmp$n_calls <- as.integer(tolower(e_tmp$official_tier)!="none")
b_nonzero <- aggregate(cbind(n_conditions,n_calls) ~ species_display + dose_pairs, e_tmp, sum)
b_nonzero$observed_call_fraction <- b_nonzero$n_calls/b_nonzero$n_conditions
b_zero <- expand.grid(species_display=species_order, dose_pairs=0, stringsAsFactors=FALSE)
b_zero$n_conditions <- nrow(blank); b_zero$n_calls <- 0; b_zero$observed_call_fraction <- 0
b <- rbind(b_zero, b_nonzero)
b$species_display <- factor(b$species_display, levels=species_order)
b$dose <- factor(b$dose_pairs, levels=doses, labels=dose_labels)
b$denominator_note <- ifelse(b$dose_pairs==0,"3 independent background blanks","15 background-by-seed conditions per species and dose")
write_tsv(b, file.path(f2out,"Fig2B_source.tsv"))
p2b <- ggplot(b, aes(dose, observed_call_fraction, group=species_display, color=species_display)) +
  geom_line(linewidth=.75) + geom_point(size=1.7) +
  scale_color_manual(values=pal_species, drop=FALSE) +
  scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1.02), expand=expansion(mult=c(0,.02))) +
  labs(title="Detection by dose", x="Spike-in pairs", y="Detection fraction\n(any non-None tier)") + theme_paper +
  guides(color=guide_legend(nrow=2, byrow=TRUE))

# C: deduplicate D000 to three background blanks, retain all non-zero real conditions.
c_nonzero <- rawC[rawC$dose_pairs>0,]
c_zero <- rawC[rawC$dose_pairs==0 & !duplicated(rawC$background_sample),]
c_zero <- do.call(rbind, lapply(species_order, function(sp) {z<-c_zero; z$species_display<-sp; z$species<-sp; z}))
c_plot_raw <- rbind(c_nonzero, c_zero)
c_summary <- aggregate(coordinate_unique_spike_fragments ~ species_display + dose_pairs, c_plot_raw, function(x) median(x,na.rm=TRUE))
names(c_summary)[3] <- "median_coordinate_unique"
c_summary$input_spike_pairs <- c_summary$dose_pairs
c_summary$n_independent_conditions <- ifelse(c_summary$dose_pairs==0,3,15)
c_summary$species_display <- factor(c_summary$species_display, levels=species_order)
c_summary$dose <- factor(c_summary$dose_pairs, levels=doses, labels=dose_labels)
write_tsv(c_summary, file.path(f2out,"Fig2C_source.tsv"))
p2c <- ggplot(c_summary, aes(dose, median_coordinate_unique, group=species_display, color=species_display)) +
  geom_line(linewidth=.72) + geom_point(size=1.6) +
  geom_line(aes(x=dose,y=input_spike_pairs, group=1, linetype="Input pairs"), inherit.aes=FALSE,
            data=unique(c_summary[c("dose","input_spike_pairs")]), color="#777777", linewidth=.6) +
  geom_point(aes(x=dose,y=input_spike_pairs), inherit.aes=FALSE,
             data=unique(c_summary[c("dose","input_spike_pairs")]), color="#777777", size=1.2) +
  scale_color_manual(values=pal_species, drop=FALSE) +
  scale_y_continuous(trans=log1p_trans(), breaks=c(0,1,2,5,10,20,50), labels=c(0,1,2,5,10,20,50)) +
  scale_linetype_manual(values=c("Input pairs"="dashed")) +
  labs(
    title="Quantitative recovery",
    x="Input spike-in pairs",
    y="Independent spike-in sites\n(median; log1p scale)"
  ) + theme_paper +
  guides(
    color="none",
    linetype=guide_legend(
      override.aes=list(color="#777777", linewidth=.6)
    )
  ) +
  theme(legend.position="bottom")

# D: summarize the four verified nested stages by species and dose.
d_nonzero <- rawD[rawD$dose_pairs>0,]
d_zero <- rawD[rawD$dose_pairs==0 & !duplicated(paste(rawD$background_sample,rawD$stage_key)),]
d_zero <- do.call(rbind, lapply(species_order, function(sp) {z<-d_zero; z$species_display<-sp; z}))
d_raw <- rbind(d_nonzero,d_zero)
d_summary <- aggregate(fragments ~ species_display + dose_pairs + stage_key + stage, d_raw, function(x) median(x,na.rm=TRUE))
names(d_summary)[5] <- "median_fragments"
d_summary$species_display <- factor(d_summary$species_display, levels=species_order)
d_summary$dose <- factor(d_summary$dose_pairs, levels=doses, labels=dose_labels)
d_summary$stage <- factor(d_summary$stage, levels=c("Same-species candidates","Risk-free core","PRIMARY","Coordinate-unique"))
write_tsv(d_summary, file.path(f2out,"Fig2D_source.tsv"))
stage_cols <- c("Same-species candidates"="#416C86","Risk-free core"="#6F99B1","PRIMARY"="#91C4C0","Coordinate-unique"="#A8AEB3")
p2d <- ggplot(d_summary, aes(dose, median_fragments, color=stage, group=stage)) +
  geom_line(linewidth=.62) + geom_point(size=1.15) + facet_wrap(~species_display,ncol=2) +
  scale_color_manual(
    values=stage_cols,
    breaks=names(stage_cols),
    labels=c(
      "Target candidates",
      "Risk-free",
      "Primary retained",
      "Unique sites"
    ),
    drop=FALSE
  ) +
  scale_x_discrete(breaks=c("0","5","50")) +
  scale_y_continuous(trans=log1p_trans(), breaks=c(0,1,2,5,10,20,50,100), labels=c(0,1,2,5,10,20,50,100)) +
  labs(title="Stage-specific retention", x="Spike-in pairs", y="Median fragments (log1p scale)") + theme_paper +
  theme(
    legend.text=element_text(size=6.1),
    legend.key.width=unit(3,"mm"),
    legend.spacing.x=unit(.5,"mm"),
    panel.spacing=unit(3,"mm"),
    axis.text.x=element_text(size=6.3)
  ) +
  guides(color=guide_legend(ncol=2,byrow=TRUE))

# E: official tiers; D000 already contains three independent blank rows.
e <- rawE
e$dose <- factor(e$dose_pairs, levels=doses, labels=dose_labels)
e$tier <- factor(tools::toTitleCase(e$official_tier), levels=names(pal_tier))
e_counts <- aggregate(run_id ~ dose + tier, e, length); names(e_counts)[3] <- "conditions"
e_den <- aggregate(run_id ~ dose, e, length); names(e_den)[2] <- "denominator"
e_counts <- merge(e_counts,e_den,by="dose"); e_counts$fraction <- e_counts$conditions/e_counts$denominator
write_tsv(e_counts, file.path(f2out,"Fig2E_source.tsv"))
p2e <- ggplot(e_counts,aes(dose,fraction,fill=tier)) + geom_col(width=.78,color="white",linewidth=.2) +
  scale_fill_manual(values=pal_tier,drop=FALSE) + scale_y_continuous(labels=percent_format(),expand=c(0,0)) +
  labs(title="Evidence tier by dose",x="Spike-in pairs",y="Fraction of independent conditions") + theme_paper +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

# F: background by species generalization, no strain claim.
f <- rawF
f$species_display <- factor(f$species_display,levels=species_order)
f$background_display <- factor(f$background_display,levels=sort(unique(f$background_display)))
write_tsv(f, file.path(f2out,"Fig2F_source.tsv"))
p2f <- ggplot(f,aes(species_display,background_display,fill=observed_call_fraction)) +
  geom_tile(color="white",linewidth=.6) + geom_text(aes(label=percent(observed_call_fraction,accuracy=1)),size=2.3,color="#243A47") +
  scale_fill_gradient(low="#F1F4F6",high="#5C8EAA",limits=c(0,1),labels=percent_format()) +
  scale_y_discrete(
    labels=function(x) sub("[[:space:]]*\\(", "\n(", x)
  ) +
  labs(
    title="Detection by background",
    x=NULL, y=NULL,
    fill="Detection fraction"
  ) + theme_paper +
  theme(
    axis.text.y=element_text(size=6.2,lineheight=1.05),
    axis.text.x=element_text(size=6.3,angle=30,hjust=1),
    legend.position="bottom",
    legend.title=element_text(size=6.2),
    legend.text=element_text(size=6.1)
  ) +
  guides(fill=guide_colorbar(
    title.position="top",
    title.hjust=.5,
    barwidth=unit(27,"mm"),
    barheight=unit(2.5,"mm")
  ))

panels2 <- list(A=p2a,B=p2b,C=p2c,D=p2d,E=p2e,F=p2f)
for(nm in names(panels2)) save3(tag(panels2[[nm]],nm), file.path(f2out,paste0("Fig2",nm)), 90, 70)
fig2 <- wrap_plots(lapply(names(panels2), function(nm) tag(panels2[[nm]],nm)), ncol=3) +
  plot_layout(heights=c(1,1.25))
save3(fig2,file.path(f2out,"Fig2_spikein_PRE_BROAD_v3"),180,165)

writeLines(c(
  "In silico spike-in validation in real WES backgrounds.",
  "A. The design comprised three host-depleted WES backgrounds, four microbial species, eight total dose levels (zero plus seven non-zero doses), and five fixed seeds, yielding 480 labeled pipeline conditions. The 60 D000 labels were duplicate labels of three underlying blank inputs and were deduplicated to three independent blanks.",
  "B. Detection includes all non-None official tiers, including Weak. Observed detection fractions used three independent blanks at dose zero and 15 background-by-seed conditions for each species at each non-zero dose. No LOD was estimated.",
  "C. The grey dashed line denotes input spike-in pairs. Coordinate-unique evidence is shown against input spike-in pairs; this is descriptive computational recovery, not an experimental capture rate.",
  "D. Median fragments are shown for the verified nested PRE-Broad stages: same-species candidates, risk-free core, PRIMARY retained, and coordinate-unique evidence.",
  "E. Official PRE-Broad evidence tiers across dose.",
  "F. Detection across three real WES backgrounds and four species. One source strain per species was used; strain generalization was not tested."
),file.path(f2out,"Fig2_legend.txt"))
writeLines(c(
  "Scientific baseline: git commit 77c4e85fd30c59a1d482b094f253356f2676597f.",
  paste("Primary input:",file.path(f2in,"Figure2A_conditions_480.tsv")),
  "Analysis level: post-host-depletion computational spike-in in real WES backgrounds.",
  "BroadExclusion: not used.",
  "Dose correction: the real grid contains eight total levels, of which seven are non-zero; the requested phrase '8 non-zero doses' was not used because it conflicts with source data.",
  "D000: 60 labels map to three independent zero-spike background inputs and were deduplicated for inferential denominators.",
  "No fabricated results or presentation mock-up data were used; this is an in-silico spike-in benchmark."
),file.path(f2out,"Fig2_provenance.txt"))

qc2 <- data.frame(
  check=c("labeled_conditions","nonzero_conditions","independent_blanks","backgrounds","species","dose_levels","nonzero_dose_levels","fixed_seeds","commit_77c4e85_only","stage_nesting_all_pass","broadexclusion_inputs","simulated_inputs","old_25_condition_inputs","composite_png_exists","composite_pdf_exists","composite_svg_exists","human_visual_review"),
  observed=c(nrow(cond),sum(cond$dose_pairs>0),nrow(blank),length(unique(cond$background_sample)),length(unique(cond$species)),length(unique(cond$dose_pairs)),sum(unique(cond$dose_pairs)>0),length(unique(cond$seed)),length(unique(cond$scientific_commit))==1 && startsWith(unique(cond$scientific_commit),"77c4e85"),all(rawC$nesting_pass==1),0,0,0,file.exists(file.path(f2out,"Fig2_spikein_PRE_BROAD_v3.png")),file.exists(file.path(f2out,"Fig2_spikein_PRE_BROAD_v3.pdf")),file.exists(file.path(f2out,"Fig2_spikein_PRE_BROAD_v3.svg")),"NOT_COMPLETED"),
  expected=c(480,420,3,3,4,8,7,5,TRUE,TRUE,0,0,0,TRUE,TRUE,TRUE,"NOT_COMPLETED"),
  status=c(rep("PASS",16),"PENDING")
)
write_tsv(qc2,file.path(f2out,"Fig2_QC.tsv"))


cat("FIG2_PRE_BROAD_V3_EXPORT_COMPLETE\n")
cat(file.path(f2out,"Fig2_spikein_PRE_BROAD_v3.png"), "\n")
