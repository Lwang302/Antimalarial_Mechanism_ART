#!/usr/bin/env Rscript



# For all users
# devtools::install_github("Tong-Chen/ImageGP")
# For chinese users
# devtools::install_git("https://gitee.com/ct5869/ImageGP.git")
library(ImageGP)
library(stringr)
library(htmlwidgets)
library(plotly)
data = "/home/Cloud_Platform/Cloud_Platform/public/system_file/enrichment/goeast.enrich.txt"
xvariable = "log_odds_ratio"
yvariable = "Term"
size_variable = "q"
color_variable = "p"
xvariable_order = "NULL"
yvariable_order = "NULL"
shape_variable_order = "NULL"
title = "Enrichment plot of GOEAST output results"
x_label = "NULL"
y_label = "NULL"
yvariable_width = 60
manual_color_vector = "Pastel1"
log10_transform_variable = "NULL"
sqrt_transform_variable = "NULL"
legend.position = "top"
xtics_angle = 0
shape_variable = "Ontology"
coordinate_flip = FALSE
scale_size_min = NULL
scale_size_max = NULL
width = 25
height = 15
outputprefix = "/home/Cloud_Platform/Cloud_Platform/public/user/_visitors/result_output/Enrichment_plot102024_01_17_11_44_30_657065/e201e143-73cf-4f26-9da6-019be6c133a4"
outputpictype = "pdf"
saveppt = FALSE
savehtml = FALSE

if (data == "") {
  script = sub(".*=", "", commandArgs()[4])
  #print(script)
  system(paste(script, "-h"))
  stop("At least -f is required!")
}



if (outputprefix == "") {
  outputprefix = data
}

filename = paste0(outputprefix, '.enrichment.', outputpictype)

xvariable_order = sp_string2vector(xvariable_order)
yvariable_order = sp_string2vector(yvariable_order)
shape_variable_order = sp_string2vector(shape_variable_order)
manual_color_vector = sp_string2vector(manual_color_vector)


cat(sp_current_time(), "Starting...\n")

sp_enrichment(
  data = data,
  xvariable = xvariable,
  yvariable = yvariable,
  size_variable = size_variable,
  color_variable = color_variable,
  xvariable_order = xvariable_order,
  yvariable_order = yvariable_order,
  shape_variable_order = shape_variable_order,
  title = title,
  x_label = x_label,
  y_label = y_label,
  yvariable_width = yvariable_width,
  manual_color_vector = manual_color_vector,
  log10_transform_variable = log10_transform_variable,
  sqrt_transform_variable = sqrt_transform_variable,
  legend.position = legend.position,
  xtics_angle = xtics_angle,
  scale_size_min = scale_size_min,
scale_size_max = scale_size_max,
  shape_variable = shape_variable,
  coordinate_flip = coordinate_flip,
  filename = filename,
  width = width,
  height = height,
  saveppt = saveppt,
  savehtml = savehtml
)
cat(sp_current_time(), "Success.\n")

