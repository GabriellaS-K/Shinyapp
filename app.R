# Load Packages
library(shiny)
library(dplyr)
library(ggplot2)
library(stringr)
library(lme4)
library(lmerTest)
library(here)
#here::here()


# Read data from RDS files

biochem_data <- readRDS(here("data/biochem.rds"))
olink_data <- readRDS(here("data/olink.rds"))
hep_data<- readRDS(here("data/hep.rds"))

biochem_data$timepoint<-as.factor(biochem_data$timepoint)
olink_data$timepoint<-as.factor(olink_data$timepoint)
hep_data$timepoint<-as.factor(hep_data$timepoint)

markers_biochem <- unique(gsub("_biochem$", "", names(biochem_data)[grepl("_biochem$", names(biochem_data))]))
markers_olink_log10 <- unique(gsub("_log10$", "", names(olink_data)[grepl("_log10$", names(olink_data))]))
markers_hep <- unique(gsub("_hep$", "", names(hep_data)[grepl("_hep$", names(hep_data))]))

######################################### UI #############################################
ui <- fluidPage(
  titlePanel("Marker Boxplots"), #Sets the title of the Shiny app

  
  sidebarLayout(
    sidebarPanel( #Sidebar inputs
      selectInput(
        inputId = "olink_select",  # ID of the first selectInput
        label = "Select a marker:",  # Label for the first selectInput
        choices = markers_olink_log10,  # Choices for the first selectInput
        selected = NULL  # Initial selected value for the first selectInput
      ),
      selectInput(
        inputId = "biochem_select",  # ID of the second selectInput
        label = "Select a biochem column:",  # Label for the second selectInput
        choices = markers_biochem,  # Choices for the second selectInput
        selected = NULL  # Initial selected value for the second selectInput
      ),
      selectInput(
        inputId = "hepcidin_select",  # ID of the third selectInput
        label = "Select hepcidin column:",  # Label for the third selectInput
        choices = markers_hep,  # Choices for the third selectInput
        selected = NULL  # Initial selected value for the third selectInput
      ),
      width = 3 #Width of the sidebar panel
    ),
    
    mainPanel(  # Contains the main panel outputs
      plotOutput(outputId = "olink_boxplot"),  # Output for the Olink boxplot plot
      textOutput(outputId = "olink_anova_result"),  # Output for the Olink ANOVA result
      plotOutput(outputId = "biochem_boxplot"),  # Output for the biochem boxplot plot
      textOutput(outputId = "biochem_anova_result"),  # Output for the biochem ANOVA result
      plotOutput(outputId = "hep_boxplot"),  # Output for the hepcidin boxplot plot
      textOutput(outputId = "hep_anova_result")  # Output for the hepcidin ANOVA result
    )
  )
)
######################################### SERVER #############################################
server <- function(input, output) {
  
  output$olink_boxplot <- renderPlot({
    selected_marker <- input$olink_select  # Get the selected marker from the first selectInput
    column_name <- paste0(selected_marker, "_log10")  # Create the column name by appending "_log10" to the selected marker
    
    if (!column_name %in% names(olink_data)) {  # Check if the column exists in the olink_data
      message("Column '", column_name, "' doesn't exist in olink data.")
      return(NULL)
    }
    
    olink_plot_data <- olink_data %>%
      select(timepoint, ID, !!sym(column_name)) %>%
      filter(!is.na(.[[column_name]]) & is.finite(.[[column_name]]))  # Filter the data, removing rows with NA or infinite values
    
    # Perform linear mixed-effects regression and calculate ANOVA
    lmer_model <- lmer(formula = as.formula(paste(column_name, "~ timepoint + (1 | ID)")), data = olink_plot_data)
    anova_result <- anova(lmer_model, type = "III")
    p_value <- formatC(anova_result[["Pr(>F)"]][1], format = "f", digits = 3)  # Format the p-value
    p_value_display <- ifelse(as.numeric(p_value) < 0.05, paste0(p_value, "*"), p_value)  # Add asterisk if p-value is significant
    
    # Create the boxplot object using ggplot
    olink_plot_object <- ggplot(data = olink_plot_data, aes(x = timepoint, y = .data[[column_name]], color = timepoint)) +
      geom_boxplot() +
      labs(x = "timepoint", y = selected_marker) +
      theme_classic()
    
    print(olink_plot_object)  # Print the boxplot object
  })
  
  output$olink_anova_result <- renderText({
    selected_marker <- input$olink_select  # Get the selected marker from the first selectInput
    column_name <- paste0(selected_marker, "_log10")  # Create the column name by appending "_log10" to the selected marker
    
    if (!column_name %in% names(olink_data)) {  # Check if the column exists in the olink_data
      message("Column '", column_name, "' doesn't exist in olink.")
      return(NULL)
    }
    
    olink_plot_data <- olink_data %>%
      select(timepoint, ID, !!sym(column_name)) %>%
      filter(!is.na(.[[column_name]]) & is.finite(.[[column_name]]))  # Filter the data, removing rows with NA or infinite values
    
    # Perform linear mixed-effects regression and calculate ANOVA
    lmer_model <- lmer(formula = as.formula(paste(column_name, "~ timepoint + (1 | ID)")), data = olink_plot_data)
    anova_result <- anova(lmer_model, type = "III")
    p_value <- formatC(anova_result[["Pr(>F)"]][1], format = "f", digits = 5)  # Format the p-value
    p_value_display <- ifelse(as.numeric(p_value) < 0.05, paste0(p_value, "*"), p_value)  # Add asterisk if p-value is significant
    
    paste("Repeated Measures ANOVA p-value:", p_value_display)  # Return the ANOVA p-value as text
  })
  
  output$biochem_boxplot <- renderPlot({
    selected_column <- paste0(input$biochem_select, "_biochem")  # Get the selected column from the second selectInput
    
    if (!selected_column %in% colnames(biochem_data)) {  # Check if the column exists in the biochem_data
      message("Column '", selected_column, "' doesn't exist in the Biochem Data.")
      return(NULL)
    }
    
    biochem_plot_data <- biochem_data %>%
      select(timepoint, ID, !!sym(selected_column)) %>%
      filter(!is.na(.[[selected_column]]) & is.finite(.[[selected_column]]))  # Filter the data, removing rows with NA or infinite values
    
    # Perform linear mixed-effects regression and calculate ANOVA
    lmer_model <- lmer(formula = as.formula(paste(selected_column, "~ timepoint + (1 | ID)")), data = biochem_plot_data)
    anova_result <- anova(lmer_model, type = "III")
    p_value <- formatC(anova_result[["Pr(>F)"]][1], format = "f", digits = 3)  # Format the p-value
    p_value_display <- ifelse(as.numeric(p_value) < 0.05, paste0(p_value, "*"), p_value)  # Add asterisk if p-value is significant
    
    # Create the boxplot object using ggplot
    biochem_plot_object <- ggplot(data = biochem_plot_data, aes(x = timepoint, y = .data[[selected_column]], color = timepoint)) +
      geom_boxplot() +
      labs(x = "timepoint", y = selected_column) +
      theme_classic()
    
    print(biochem_plot_object)  # Print the boxplot object
  })
  
  output$biochem_anova_result <- renderText({
    selected_column <- paste0(input$biochem_select, "_biochem")  # Get the selected column from the second selectInput
    
    if (!selected_column %in% colnames(biochem_data)) {  # Check if the column exists in the biochem_data
      message("Column '", selected_column, "' doesn't exist in the Biochem Data.")
      return(NULL)
    }
    
    biochem_plot_data <- biochem_data %>%
      select(timepoint, ID, !!sym(selected_column)) %>%
      filter(!is.na(.[[selected_column]]) & is.finite(.[[selected_column]]))  # Filter the data, removing rows with NA or infinite values
    
    # Perform linear mixed-effects regression and calculate ANOVA
    lmer_model <- lmer(formula = as.formula(paste(selected_column, "~ timepoint + (1 | ID)")), data = biochem_plot_data)
    anova_result <- anova(lmer_model, type = "III")
    p_value <- formatC(anova_result[["Pr(>F)"]][1], format = "f", digits = 5)  # Format the p-value
    p_value_display <- ifelse(as.numeric(p_value) < 0.05, paste0(p_value, "*"), p_value)  # Add asterisk if p-value is significant
    
    paste("Repeated Measures ANOVA p-value for Biochem Data:", p_value_display)  # Return the ANOVA p-value as text
  })
  
  output$hep_boxplot <- renderPlot({
    selected_marker <- input$hepcidin_select  # Get the selected marker from the third selectInput
    column_name <- paste0(selected_marker, "_hep")  # Create the column name by appending "_hep" to the selected marker
    
    if (!column_name %in% names(hep_data)) {  # Check if the column exists in the hep_data
      message("Column '", column_name, "' doesn't exist in hep data.")
      return(NULL)
    }
    
    hep_plot_data <- hep_data %>%
      select(timepoint, ID, !!sym(column_name)) %>%
      filter(!is.na(.[[column_name]]) & is.finite(.[[column_name]]))  # Filter the data, removing rows with NA or infinite values
    
    # Perform linear mixed-effects regression and calculate ANOVA
    lmer_model <- lmer(formula = as.formula(paste(column_name, "~ timepoint + (1 | ID)")), data = hep_plot_data)
    anova_result <- anova(lmer_model, type = "III")
    p_value <- formatC(anova_result[["Pr(>F)"]][1], format = "f", digits = 3)  # Format the p-value
    p_value_display <- ifelse(as.numeric(p_value) < 0.05, paste0(p_value, "*"), p_value)  # Add asterisk if p-value is significant
    
    # Create the boxplot object using ggplot
    hep_plot_object <- ggplot(data = hep_plot_data, aes(x = timepoint, y = .data[[column_name]], color = timepoint)) +
      geom_boxplot() +
      labs(x = "timepoint", y = selected_marker) +
      theme_classic()
    
    print(hep_plot_object)  # Print the boxplot object
  })
  
  output$hep_anova_result <- renderText({
    selected_marker <- input$hepcidin_select  # Get the selected marker from the third selectInput
    column_name <- paste0(selected_marker, "_hep")  # Create the column name by appending "_hep" to the selected marker
    
    if (!column_name %in% names(hep_data)) {  # Check if the column exists in the hep_data
      message("Column '", column_name, "' doesn't exist in hep.")
      return(NULL)
    }
    
    hep_plot_data <- hep_data %>%
      select(timepoint, ID, !!sym(column_name)) %>%
      filter(!is.na(.[[column_name]]) & is.finite(.[[column_name]]))  # Filter the data, removing rows with NA or infinite values
    
    # Perform linear mixed-effects regression and calculate ANOVA
    lmer_model <- lmer(formula = as.formula(paste(column_name, "~ timepoint + (1 | ID)")), data = hep_plot_data)
    anova_result <- anova(lmer_model, type = "III")
    p_value <- formatC(anova_result[["Pr(>F)"]][1], format = "f", digits = 5)  # Format the p-value
    p_value_display <- ifelse(as.numeric(p_value) < 0.05, paste0(p_value, "*"), p_value)  # Add asterisk if p-value is significant
    
    paste("Repeated Measures ANOVA p-value:", p_value_display)  # Return the ANOVA p-value as text
  })
  
}


shinyApp(ui = ui, server = server)



