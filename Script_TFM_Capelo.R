library(tidyverse)
library(gt)
library(gtsummary)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(car)
library(caret)
library(pROC)
library(glmnet)
library(gbm)

# Limpieza de la base de datos original

setwd('~/UNIR/TFM/TFM_Capelo') #Establecer aquí el directorio clonado
db_original= read.csv("base_datos_original(MSK, Nat Med 2024).csv", sep= '\t') #Establecer aquí la ubicación del archivo donde se guarda la base datos original

## Renombramos las varibles sustituyendo los punto por _ y eliminando paréntesis
colnames(db_original) <- gsub('\\.', '_', colnames(db_original)) #Sustituir espacios por _ en nombre de variables
colnames(db_original) <- gsub ('\\(|\\)', '', colnames(db_original)) #Eliminar paréntesis en nombre de variables
db_original <- db_original %>%
  rename( Overall_Survival_Months = Overall_Survival__Months_, Diagnosis_to_Collection_Year= Diagnosis_to_Collection__Year_) %>%
  mutate(Age_at_Diagnosis = gsub('>','', Age_at_Diagnosis))

## Asignamos la clase correspondiente a cada una de las variables
var_numericas <- c('Age_at_Diagnosis','Diagnosis_to_Collection_Year','Fraction_Genome_Altered','KRAS_Mutation_Dosage','KRAS_Mutant_Copies_Gained','Mutation_Count','Overall_Survival_Months','Number_of_Samples_Per_Patient','Sample_coverage')
var_logical <- c('Evaluable_for_Allelic_Imbalance','Ashkenazi_Jewish_Ancestry','Autoimmune_Disease_History','Coronary_Artery_Disease_History','Cancer_History','Diabetes_History','Hypertension_Hisotry','Pancreatitis_History','Resection_Surgery','Systemic_Treatment_1','Targeted_Therapy')
var_character <- c('Study_ID','Patient_ID','Sample_ID')
var_factor <- c('Adjuvant_Treatment', 'KRAS_Allele_Selection','Cancer_Type', 'Cancer_Type_Detailed','Copy_Number_State', 'Disease_Status_at_Sequencing', 'Gene_Panel', 'Genomic_Subtype', 'Sub_threshold_KRAS_Mutation', 'KRAS_Somatic_Variant', 'MSI_Type', 'Neoadjuvant_Treatment', 'Oncotree_Code', 'Overall_Survival', 'Sample_Type', 'KRAS_Allele_Selection_State', 'Sex', 'Somatic_Status',  'Systemic_Treatment', 'Systemic_Treatment_1','Systemic_Treatment_2', 'Systemic_Treatment_3', 'Systemic_Treatment_1', 'Systemic_Treatment_4', 'Systemic_Treatment_5', 'TMB_H', 'Smoking_History', 'Tumor_Location', 'Whole_Genome_Doubling_Status')
db <- db_original%>%
  filter(Overall_Survival == '1:DECEASED')
for (x in var_numericas){db[[x]] <- as.numeric(db[[x]])}
for (x in var_logical){db[[x]] <- as.logical(db[[x]])}
for (x in var_character){db[[x]] <- as.character(db[[x]])}
for (x in var_factor){db[[x]] <- as.factor(db[[x]])}

## Abreviamos determinadas categorías para hacer más manejables las tablas y gráficas.
db <- db %>%
  mutate(Institute_Source= case_when(Institute_Source == 'Memorial Sloan Kettering Cancer Center' ~ 'MSKCC',
                                     Institute_Source == 'Mount Sinai St. Lukes' ~ 'MSSL',
                                     Institute_Source == 'Saint Francis Hospital' ~ 'SFH',
                                     Institute_Source == 'Kings County Hospital Center,Pathology Department,451 Clarkson Avenue,Brooklyn Ny 11203,718 245 5374' ~ 'KCHC',
                                     Institute_Source == 'New York Hospital Queens' ~ 'NYHQ',
                                     Institute_Source == 'Mount Sinai Laboratory' ~ 'MSL',
                                     Institute_Source == 'The Mount Sinai Medical Center' ~ 'MSC',
                                     Institute_Source == 'Other Institute' ~ 'Other')) %>%
  mutate(Stage_at_Diagnosis= case_when(
    Stage_at_Diagnosis == 'Borderline Resectable/Locally Advanced' ~ 'BR/LA',
    TRUE ~ Stage_at_Diagnosis)) %>% 
  mutate(Ancestry= case_when(Ancestry == 'Admixed/other' ~ 'Other',
                             Ancestry == 'African (AFR)' ~ 'AFR',
                             Ancestry == 'East Asian (EAS)' ~ 'EAS',
                             Ancestry == 'European (EUR)' ~ 'EUR',
                             Ancestry == 'South Asian (SAS)' ~ 'SAS'
  ))

## Generamos dataframe informativo de variables (Incluyendo clase y %NA)
var_df <- data.frame()
for (variable in names(db)){
  na_porc <- round(mean(is.na(db[[variable]]))*100, 2)
  tipo <- class(db[[variable]])
  var_df <- rbind(var_df, list(gsub('_', ' ',variable),  tipo, na_porc))
}
colnames(var_df) <- c('Variable', 'Tipo', 'NA (%)')

## Importamos descripción de las variables y las añadimos al dataframe
descripciones <- c()
description_var <- readLines('descripción_variables.txt')
for (y in description_var){
  d <- strsplit(y, ':')
  descripciones <- c(descripciones, trimws(d[[1]][2])) #trimws para eliminar espacios delante del string
}
var_df <- cbind(var_df, 'Descripción' = as.matrix(descripciones))

## Depuración base de datos
### Eliminamos variables con >10% de NA, y aquellas entradas incompletas
variables_eliminar <- var_df %>%
  filter(`NA (%)` > 10)
db <- db %>%
  select(- gsub(' ', '_', variables_eliminar$Variable))
num_pacientes_na <- sum(apply(is.na(db), 1, any))
db <- na.omit(db) #elimanmos aquellas entradas que no estén completas

### Gestión de outliers 
patients_outliers <- c()
for (variable in names(db)){
  if(class(db[[variable]]) %in% c('integer', 'numeric')){
    iqr= IQR(db[[variable]])
    q1 = quantile(db[[variable]], 0.25)
    q3 = quantile(db[[variable]], 0.75)
    lower_limit = q1 - 1.5 * iqr
    upper_limit = q3 + 1.5 * iqr
    outliers <- db %>%
      filter(db[[variable]] < lower_limit | db[[variable]] > upper_limit)
    patients_outliers <- unique(c(patients_outliers, outliers$Patient_ID))}}
length(patients_outliers)
db <- db %>%
  filter(!(db$Patient_ID %in% patients_outliers))

### Calculamos la varianza y número de categorías para las variables restantes
varianzas <- c()
categorias <- c()
for (variable in var_df$Variable){
  if(!(variable %in% variables_eliminar$Variable)){
    variable <- gsub(' ','_',variable)
    if (class(db[[variable]]) == 'numeric'){
      varianza <- round(var(db[[variable]], na.rm=TRUE), 2)
      categories <- '-'  }
    else { varianza <- '-'
    categories <- length(unique(db[[variable]]))}}
  else{
    varianza <- '-'
    categories <- '-'}
  varianzas <- c(varianzas, varianza)
  categorias <- c(categorias, categories)
}
var_df <- cbind(var_df, 'Varianza'= varianzas, 'Categorías'= categorias)

### Eliminamos aquellas entradas con varianza 0 y solo 1 categoría del dataframe db 
variables_eliminar_2 <- var_df %>%
  filter(Varianza == 0 | Categorías == 1)
db <- db %>%
  select(- gsub(' ', '_', variables_eliminar_2$Variable))

## Convetimos en tabla y exportamos el dataframe con las descripciones de las variables
var_tabla <- var_df %>%
  gt() %>%
  tab_spanner(
    label= 'Tras eliminación NA',
    columns = c(Varianza, Categorías)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners(spanners = "Tras eliminación NA")) %>%
  tab_style(
    style = cell_text(weight = "bold"),  # Aplica negrita
    locations = cells_column_labels() 
  ) %>%
  tab_style(
    style = list (cell_fill(color='#FFBEB2')),
    locations = cells_body(rows= Variable %in% c(variables_eliminar$Variable, variables_eliminar_2$Variable))
  )

##Dividimos el DataFrame en grupos
median_survival <- median(db$Overall_Survival_Months)
db <- db %>%
  mutate(Group_Survival= case_when(
    Overall_Survival_Months >= median_survival ~ 'HIGHER',
    Overall_Survival_Months < median_survival ~ 'LOWER'
  ))
db$Group_Survival <- as.factor(db$Group_Survival)




# Análisis comparativo descriptivo de la muestra
db_descriptivo <- db %>%
  select(- c(Sample_ID, Patient_ID))

## Análisis de normalidad (q-q plot)
qq_plots1 <- list()
qq_plots2 <- list()
qq_contador <- 1
for (variable in names(db_descriptivo)){
  # Si es continua
  if(class(db_descriptivo[[variable]]) %in% c('integer', 'numeric')){
    for(grupo in c('HIGHER','LOWER')){
      db_grupo <- db_descriptivo %>%
        filter(Group_Survival == grupo)
      #definimos etiquetas según que grupo sea
      if(qq_contador < 3){
        titulo <- grupo
      }
      else{titulo <- ''}
      if(grupo == 'HIGHER'){
        mi_labs <- labs(title= titulo, y= gsub('_',' ', variable), x='Cuantiles normales')}
      else if(grupo == 'LOWER'){
        mi_labs <- list(labs(title= titulo, y= 'Cuantiles observados', x='Cuantiles normales'), 
                        scale_y_continuous(position='right'))}
      qqplot <- ggplot(db_grupo, aes(sample = db_grupo[[variable]])) +
        geom_qq() +
        geom_qq_line(color = "red") +
        mi_labs +
        theme_minimal() 
      
      if(qq_contador <= 8){
        qq_plots1[[qq_contador]] <- ggplotGrob(qqplot)
        qq_contador <- qq_contador + 1}
      else if(qq_contador > 8) {qq_plots2[[qq_contador-8]] <- ggplotGrob(qqplot)
      qq_contador <- qq_contador + 1}}}}

qq_arrange1 <- grid.arrange(grobs= qq_plots1, ncol= 2)
qq_arrange2 <- grid.arrange(grobs= qq_plots2, ncol= 2)

## Tabla descriptiva
normales <- c('Age at Diagnosis', 'Mutation Count', 'Sample coverage', 'TMB  nonsynonymous ')
db_tabla <- db %>%
  select(- c('Patient_ID', 'Sample_ID'))
colnames(db_tabla) <- gsub('_', ' ', colnames(db_tabla))

no_normales <- db_tabla%>%
  select_if(is.numeric) %>%
  select(-normales)%>%
  colnames()
### Tabla de totales
total_table <- db_tabla %>%
  select(-c('Group Survival'))%>%
  tbl_summary(
    type= list(`Mutation Count` ~ 'continuous'), 
    statistic =list(
      normales ~ '{mean} ({sd})', 
      no_normales ~ '{median} ({p25} - {p75})'),
    digits= list(all_categorical() ~ c (0,1),
                 all_continuous()~ c(2,1)))%>%
  modify_header(label='**Variable**')%>%
  bold_labels()

### Tabla higher-lower
higherlower_table <- db_tabla %>%
  filter(!(`Group Survival` == 'ALIVE'))%>%
  tbl_summary(by = 'Group Survival',
              type= list(`Mutation Count` ~ 'continuous'), 
              statistic = list(normales ~ '{mean} ({sd})', 
                               no_normales ~ '{median} ({p25} - {p75})'),
              digits= list(all_categorical() ~ c (0,1),
                           all_continuous()~ c(2,1)))%>%
  add_p(test= list(
    normales ~'t.test',
    no_normales ~ 'wilcox.test',
    all_categorical()~'chisq.test'))%>%
  
  modify_header(label='**Variable**', p.value= '**p valor**')%>%
  bold_labels()

tabla_descriptiva <- tbl_merge(list(total_table,higherlower_table), tab_spanner= c('**TOTAL**', '**HIGHER-LOWER**'))

## Representaciones gráficas
num_graph1 <- list()
cat_graph1 <-list()
num_contador <- 1
cat_contador <- 1

### Generamos los gráficos

for(variab in setdiff(names(db), c('Study_ID', 'Patient_ID', 'Sample_ID','Group_Survival'))){ #Para excluir las variables de la lista
  if(class(db_descriptivo[[variab]]) %in% c('numeric', 'integer')){
    if(gsub('_', ' ', variab) %in% normales){
      metodo <- 't.test'}
    else{metodo <- 'wilcox.test'}
    if(variab == 'Overall_Survival_Months'){
      num_graph <- ggplot(db_descriptivo, aes(x= Group_Survival, y= !!sym(variab), fill= Group_Survival)) +
        geom_boxplot(width=0.6)+
        labs(title= '', y= gsub('_', ' ', variab), x= '', y='')+
        scale_fill_manual(values= c('#EFB663', '#FFFF37'))+
        theme_light()+
        theme(legend.position = "none")+
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      num_graph1[[num_contador]] <- ggplotGrob(num_graph)
      num_contador <- num_contador +1
    }
    else{
      num_graph <- ggplot(db_descriptivo, aes(x= Group_Survival, y= !!sym(variab), fill= Group_Survival)) +
        geom_boxplot(width=0.6)+
        labs(title= '', y= gsub('_', ' ', variab), x= '', y='')+
        scale_fill_manual(values= c('#EFB663', '#FFFF37'))+
        theme_light()+
        stat_compare_means(aes(group= Group_Survival), 
                           label = "p.signif", 
                           method= metodo,
                           comparisons = list(c('HIGHER','LOWER'))) +
        theme(legend.position = "none")+
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      num_graph1[[num_contador]] <- ggplotGrob(num_graph)
      num_contador <- num_contador +1
    }}
  else{
    df_grafico_variab <- db_descriptivo %>%
      group_by(Group_Survival, !!sym(variab))%>% #Agrupamos cada grupo con cada categoría
      summarise(n=n(), .groups = 'drop')%>% #De aquí obtenemos una columna n que representa cada uno de los totales. .groups=`drop'para eliminar dicha agrupación.
      group_by(Group_Survival)%>% #Ahora agrupamos la tabla por grupo para que la suma de n sea el total de cada grupo
      mutate(Porcentaje= n/sum(n)*100) #Dividimos cada n entre el total de cada grupo
    
    #Para obtener el p valor de la prueba chi cuadrado: 
    #Función para converir el valor de p valor en código de asteriscos
    p_to_stars <- function(p) {
      if (p < 0.0001) return("****")
      else if (p < 0.001) return("***")
      else if (p < 0.01) return("**")
      else if (p < 0.05) return("*")
      else return("ns")  # no significativo
    }
    p_valor <- higherlower_table$table_body %>%
      filter(label==gsub('_', ' ', variab))%>%
      pull(p.value)
    p_label <- p_to_stars(p_valor)
    pos_centro_categorías <- (length((unique(db_descriptivo[[variab]])))+1)/2
    cat_graph <- ggplot(df_grafico_variab, aes(x=!!sym(variab),y= Porcentaje, fill= Group_Survival)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.6)+
      labs(y= paste(gsub('_', ' ', variab),' (%)'), x= '')+
      scale_fill_manual(values= c('#EFB663', '#FFFF37'))+
      scale_y_continuous(limits = c(0, 100)) +
      theme_light()+
      theme(legend.position = "none")+
      theme(plot.title = element_text(hjust = 0.5, face='bold', size=7))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7))+
      annotate("text", x = pos_centro_categorías, y = 98,
               label = p_label,
               size = 4)
    cat_graph1[[cat_contador]] <- ggplotGrob(cat_graph)
    cat_contador <- cat_contador +1
  }
}

### Generamos uno de los gráficos para obtener la leyenda
leyenda_graph <- ggplot(df_grafico_variab, aes(x=Stage_at_Diagnosis,y= Porcentaje, fill= Group_Survival)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6)+
  labs(y= paste(gsub('_', ' ', variab),' (%)'), x= '', fill='Grupo')+
  scale_fill_manual(values= c('#EFB663', '#FFFF37'))+
  scale_y_continuous(limits = c(0, 100)) +
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7))+
  theme(plot.margin= margin(t = 100, r = 10, b = 10, l = 10))
cat_leyenda <- get_legend(leyenda_graph)
cat_graph1[[12]] <- cat_leyenda

graph_num1 <- grid.arrange(grobs= num_graph1, ncol= 4)
graph_cat1 <- grid.arrange(grobs= cat_graph1, ncol=4)

#Construcción de modelos predictivos

set.seed(123)
colores <- rainbow(5)
db_hl <- db %>%
  select(- Overall_Survival_Months, - Patient_ID, -Sample_ID) #Eliminamos las variables que están directamente relacionadas con la variable dependiente y el identificador
db_hl$Group_Survival <- factor(db_hl$Group_Survival, levels = c('LOWER', 'HIGHER')) #Nos aseguramos de que LOWER sea la categoría positiva 

## Regularización tipo LASSO

hlc_lasso_x <- model.matrix(Group_Survival ~., data= db_hl)[,-1] #Convertimos las variables factor/categóricas en variables dummies (codificadas en k-1 columnas, k siendo las categorías). [,-1] para eliminar intercept
hlc_lasso_y <- db_hl$Group_Survival
hlc_lasso <- cv.glmnet(hlc_lasso_x, hlc_lasso_y, family='binomial', alpha=1) #Modelo lasso
hlc_lasso_coef <- as.data.frame(as.matrix(coef(hlc_lasso, s='lambda.min'))) #Obtenemos coeficientes
hlc_lasso_variab <- hlc_lasso_coef %>% 
  filter(s1 != 0)#Seleccionamos solo aquellas cuyos coeficientes sean distintos de 0.
hlc_lasso_variab_names <- rownames(hlc_lasso_variab)[-1] #Guardamos sus nombres

hlc_lasso_data <- as.data.frame(hlc_lasso_x[,hlc_lasso_variab_names], drop=FALSE) #filtramos el dataset con variables dummies para tener solo las variables relevantes.
hlc_lasso_data <- cbind(hlc_lasso_data, db_hl$Group_Survival)#Anadimos los valores de clasificación 
colnames(hlc_lasso_data)[colnames(hlc_lasso_data) == 'db_hl$Group_Survival'] <- 'Group_Survival' #Para cambiarle el nombre

## Dividimos el dataset en training y testing
hlc_trainIndex <- createDataPartition(hlc_lasso_data$Group_Survival, p=0.8, list=FALSE)
train_hlc <- hlc_lasso_data[hlc_trainIndex, ]
test_hlc <- hlc_lasso_data[-hlc_trainIndex,]

## Creamos una seed para garantizar la reproducibilidad de los modelos
n_folds <- 10
n_tune <- 20
seeds_list <- vector(mode = "list", length = n_folds + 1) # Lista de semillas para cada resample y el modelo final
for(i in 1:n_folds){
  seeds_list[[i]] <- sample.int(1000, n_tune)
}
seeds_list[[n_folds + 1]] <- sample.int(1000, 1)# Semilla para el modelo final
hlc_cv <- trainControl(method = "cv", number = 10, seeds= seeds_list) #definimos la validación cruzada que se llevará a cabo en los modelos. 


## kNN vecinos más cercanos 
knn_hlc <- train(Group_Survival~., 
                 data= train_hlc, 
                 method= 'knn',
                 preProcess= c('center', 'scale'),
                 trControl= hlc_cv, 
                 tuneLength = 20
)

knn_hlc_pred <- predict(knn_hlc, newdata = test_hlc)
knn_hlc_eval <- confusionMatrix(knn_hlc_pred, test_hlc$Group_Survival)
knn_hlc_prob <- predict(knn_hlc, newdata = test_hlc, type='prob')

knn_hlc_precision <- knn_hlc_eval$overall['Accuracy'] #Precisión general del modelo (TP+TN/TOTAL)
knn_hlc_sensitivity <- knn_hlc_eval$byClass['Sensitivity'] #(TAMBIÉN RECALL) todos los positivos reales, cuantos detecta (TP/FN+TP)
knn_hlc_specificity <- knn_hlc_eval$byClass['Specificity'] #De todos los NEGATIVOS reales cuantos detecta (TN/TN+FP)
knn_hlc_f1 <- knn_hlc_eval$byClass['F1']
knn_hlc_roc <- roc(test_hlc$Group_Survival, knn_hlc_prob[,2], levels = c('HIGHER', 'LOWER'))
knn_hlc_auc<- auc(knn_hlc_roc)
knn_hlc_legend <- paste('AUC k-NN:', round(knn_hlc_auc,2))
knn_hlc_metrics <- list('Precisión'= round(knn_hlc_precision[[1]],3), 'Sensibilidad' = round(knn_hlc_sensitivity[[1]],3), 'Especificidad'=round(knn_hlc_specificity[[1]],3), 'F1-score'=round(knn_hlc_f1[[1]],3), 'AUC'=round(knn_hlc_auc[[1]],3))
knn_hlc_metrics_df <- data.frame(Métricas=names(knn_hlc_metrics), Valores = unlist(knn_hlc_metrics))

### Curva ROC
plot(knn_hlc_roc, col= colores[[1]], main= 'CURVA ROC (k-NN)', ylab='Sensibilidad', xlab='Especificidad', lwd=2)
legend('bottomright', legend= knn_hlc_legend, col= colores[[1]], lwd=2)


## Suport Vector Machine SVM Linear
svm_hlc <- train(Group_Survival~., 
                 data= train_hlc, 
                 method= 'svmLinear',
                 preProcess= c('center', 'scale'),
                 trControl= hlc_cv, 
                 tuneLength = 20, 
                 prob.model= TRUE
)

svm_hlc_pred <- predict(svm_hlc, newdata = test_hlc)
svm_hlc_eval <- confusionMatrix(svm_hlc_pred, test_hlc$Group_Survival)
svm_hlc_prob <- predict(svm_hlc, newdata = test_hlc, type='prob')

svm_hlc_precision <- svm_hlc_eval$overall['Accuracy'] #Precisión general del modelo (TP+TN/TOTAL) 
svm_hlc_sensitivity <- svm_hlc_eval$byClass['Sensitivity'] #(TAMBIÉN RECALL) todos los positivos reales, cuantos detecta (TP/FN+TP) 
svm_hlc_specificity <- svm_hlc_eval$byClass['Specificity'] #De todos los NEGATIVOS reales cuantos detecta (TN/TN+FP) 
svm_hlc_f1 <- svm_hlc_eval$byClass['F1'] 
svm_hlc_metrics <- list('Support vector machine lineal', round(svm_hlc_precision[[1]],3), round(svm_hlc_sensitivity[[1]],3), round(svm_hlc_specificity[[1]],3), round(svm_hlc_f1[[1]],3)) 
svm_hlc_roc <- roc(test_hlc$Group_Survival, svm_hlc_prob[,2], levels = c('HIGHER', 'LOWER'))
svm_hlc_auc<- auc(svm_hlc_roc)
svm_hlc_legend <- paste('AUC SVM-Lineal:', round(svm_hlc_auc,2))
svm_hlc_metrics <- list('Precisión'= round(svm_hlc_precision[[1]],3), 'Sensibilidad' = round(svm_hlc_sensitivity[[1]],3), 'Especificidad'=round(svm_hlc_specificity[[1]],3), 'F1-score'=round(svm_hlc_f1[[1]],3), 'AUC'=round(svm_hlc_auc[[1]],3))
svm_hlc_metrics_df <- data.frame(Métricas=names(svm_hlc_metrics), Valores = unlist(svm_hlc_metrics))

###Curva ROC
plot(svm_hlc_roc, col= colores[[2]], main= 'CURVA ROC (SVM-Lineal)', ylab='Sensibilidad', xlab='Especificidad', lwd=2)
legend('bottomright', legend= svm_hlc_legend, col= colores[[2]], lwd=2)

##Suport Vector Machine SVM Gaussian 
svmG_hlc <- train(Group_Survival~., 
                  data= train_hlc, 
                  method= 'svmRadial', 
                  preProcess= c('center', 'scale'), 
                  trControl= hlc_cv, 
                  tuneLength = 20, 
                  prob.model=TRUE
) 

svmG_hlc_pred <- predict(svmG_hlc, newdata = test_hlc) 
svmG_hlc_eval <- confusionMatrix(svmG_hlc_pred, test_hlc$Group_Survival) 
svmG_hlc_prob <- predict(svmG_hlc, newdata = test_hlc, type='prob') 

svmG_hlc_precision <- svmG_hlc_eval$overall['Accuracy'] #Precisión general del modelo (TP+TN/TOTAL)
svmG_hlc_sensitivity <- svmG_hlc_eval$byClass['Sensitivity'] #(TAMBIÉN RECALL) todos los positivos reales, cuantos detecta (TP/FN+TP) 
svmG_hlc_specificity <- svmG_hlc_eval$byClass['Specificity'] #De todos los NEGATIVOS reales cuantos detecta (TN/TN+FP) 
svmG_hlc_f1 <- svmG_hlc_eval$byClass['F1'] 
svmG_hlc_metrics <- list('Support vector machine gaussiano', round(svmG_hlc_precision[[1]],3), round(svmG_hlc_sensitivity[[1]],3), round(svmG_hlc_specificity[[1]],3), round(svmG_hlc_f1[[1]],3)) 
svmG_hlc_roc <- roc(test_hlc$Group_Survival, svmG_hlc_prob[,2], levels = c('HIGHER', 'LOWER'))
svmG_hlc_auc<- auc(svmG_hlc_roc)
svmG_hlc_legend <- paste('SVM-Gaussiano:', round(svmG_hlc_auc,2))
svmG_hlc_metrics <- list('Precisión'= round(svmG_hlc_precision[[1]],3), 'Sensibilidad' = round(svmG_hlc_sensitivity[[1]],3), 'Especificidad'=round(svmG_hlc_specificity[[1]],3), 'F1-score'=round(svmG_hlc_f1[[1]],3), 'AUC'=round(svmG_hlc_auc[[1]],3))
svmG_hlc_metrics_df <- data.frame(Métricas=names(svmG_hlc_metrics), Valores = unlist(svmG_hlc_metrics))

### Curva ROC
plot(svmG_hlc_roc, col= colores[[3]], main= 'CURVA ROC (SVM-Gaussiano)', ylab='Sensibilidad', xlab='Especificidad', lwd=2)
legend('bottomright', legend= svmG_hlc_legend, col= colores[[3]], lwd=2)

# Random Forest
rf_hlc <- train(Group_Survival~., 
                data= train_hlc, 
                method= 'rf', 
                preProcess= c('center', 'scale'), 
                trControl= hlc_cv, 
                tuneLength = 20) 

rf_hlc_pred <- predict(rf_hlc, newdata = test_hlc) 
rf_hlc_eval <- confusionMatrix(rf_hlc_pred, test_hlc$Group_Survival) 
rf_hlc_prob <- predict(rf_hlc, newdata = test_hlc, type='prob')  

rf_hlc_precision <- rf_hlc_eval$overall['Accuracy'] #Precisión general del modelo (TP+TN/TOTAL)
rf_hlc_sensitivity <- rf_hlc_eval$byClass['Sensitivity'] #(TAMBIÉN RECALL) todos los positivos reales, cuantos detecta (TP/FN+TP) 
rf_hlc_specificity <- rf_hlc_eval$byClass['Specificity'] #De todos los NEGATIVOS reales cuantos detecta (TN/TN+FP) 
rf_hlc_f1 <- rf_hlc_eval$byClass['F1'] 
rf_hlc_metrics <- list('Random forest', round(rf_hlc_precision[[1]],3), round(rf_hlc_sensitivity[[1]],3), round(rf_hlc_specificity[[1]],3), round(rf_hlc_f1[[1]],3)) 
rf_hlc_roc <- roc(test_hlc$Group_Survival, rf_hlc_prob[,2], levels = c('HIGHER', 'LOWER'))
rf_hlc_auc<- auc(rf_hlc_roc)
rf_hlc_legend <- paste('AUC RF:', round(rf_hlc_auc,2))
rf_hlc_metrics <- list('Precisión'= round(rf_hlc_precision[[1]],3), 'Sensibilidad' = round(rf_hlc_sensitivity[[1]],3), 'Especificidad'=round(rf_hlc_specificity[[1]],3), 'F1-score'=round(rf_hlc_f1[[1]],3), 'AUC'=round(rf_hlc_auc[[1]],3))
rf_hlc_metrics_df <- data.frame(Métricas=names(rf_hlc_metrics), Valores = unlist(rf_hlc_metrics))

#Curva ROC
plot(rf_hlc_roc, col= colores[[4]], main= 'CURVA ROC (RF)', ylab='Sensibilidad', xlab='Especificidad', lwd=2)
legend('bottomright', legend= rf_hlc_legend, col= colores[[4]], lwd=2)

# Gradient Boosting Machine
gbm_hlc <- train(Group_Survival~., 
                 data= train_hlc, 
                 method= 'gbm', 
                 preProcess= c('center', 'scale'), 
                 trControl= hlc_cv, 
                 tuneLength = 20
)  

gbm_hlc_pred <- predict(gbm_hlc, newdata = test_hlc) 
gbm_hlc_eval <- confusionMatrix(gbm_hlc_pred, test_hlc$Group_Survival) 
gbm_hlc_prob <- predict(gbm_hlc, newdata = test_hlc, type='prob')  

gbm_hlc_precision <- gbm_hlc_eval$overall['Accuracy'] #Precisión general del modelo (TP+TN/TOTAL)
gbm_hlc_sensitivity <- gbm_hlc_eval$byClass['Sensitivity'] #(TAMBIÉN RECALL) todos los positivos reales, cuantos detecta (TP/FN+TP) 
gbm_hlc_specificity <- gbm_hlc_eval$byClass['Specificity'] #De todos los NEGATIVOS reales cuantos detecta (TN/TN+FP) 
gbm_hlc_f1 <- gbm_hlc_eval$byClass['F1'] 
gbm_hlc_metrics <- list('Gradient boosting machine', round(gbm_hlc_precision[[1]],3), round(gbm_hlc_sensitivity[[1]],3), round(gbm_hlc_specificity[[1]],3), round(gbm_hlc_f1[[1]],3)) 
gbm_hlc_roc <- roc(test_hlc$Group_Survival, gbm_hlc_prob[,2], levels = c('HIGHER', 'LOWER'))
gbm_hlc_auc<- auc(gbm_hlc_roc)
gbm_hlc_legend <- paste('AUC GBM:', round(gbm_hlc_auc,2))
gbm_hlc_metrics <- list('Precisión'= round(gbm_hlc_precision[[1]],3), 'Sensibilidad' = round(gbm_hlc_sensitivity[[1]],3), 'Especificidad'=round(gbm_hlc_specificity[[1]],3), 'F1-score'=round(gbm_hlc_f1[[1]],3), 'AUC'=round(gbm_hlc_auc[[1]],3))
gbm_hlc_metrics_df <- data.frame(Métricas=names(gbm_hlc_metrics), Valores = unlist(gbm_hlc_metrics))

### Curva ROC
plot(gbm_hlc_roc, col= colores[[5]], main= 'CURVA ROC (GBM)', ylab='Sensibilidad', xlab='Especificidad', lwd=2)
legend('bottomright', legend= gbm_hlc_legend, col= colores[[5]], lwd=2)

## TOP10 Variables relevantes
gbm_topvar <- varImp(gbm_hlc)
gbm_top10 <- gbm_topvar$importance %>%
  arrange(desc(Overall)) %>%
  slice(1:5)%>%
  rename('Relevancia' = Overall)%>% 
  mutate(Pos=(1:10))%>% 
  rownames_to_column(var = "Variables") %>%
  select(Pos, everything())
