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
stfr_data<- readRDS(here("data/stfr.rds"))

# Make sure timepoint is a factor
biochem_data$timepoint<-as.factor(biochem_data$timepoint)
olink_data$timepoint<-as.factor(olink_data$timepoint)
hep_data$timepoint<-as.factor(hep_data$timepoint)
stfr_data$timepoint<-as.factor(stfr_data$timepoint)

# Extracting unique marker names from column names in biochem_data that end with "_biochem"
markers_biochem <- unique(gsub("_biochem$", "", names(biochem_data)[grepl("_biochem$", names(biochem_data))]))

# Extracting unique marker names from column names in olink_data that end with "_log10"
markers_olink_log10 <- unique(gsub("_log10$", "", names(olink_data)[grepl("_log10$", names(olink_data))]))

# Extracting unique marker names from column names in hep_data that end with "_hep"
markers_hep <- unique(gsub("_hep$", "", names(hep_data)[grepl("_hep$", names(hep_data))]))

# Extracting unique marker names from column names in stfr_data that end with "_stfr"
markers_stfr <- unique(gsub("_stfr$", "", names(stfr_data)[grepl("_stfr$", names(stfr_data))]))


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
      selectInput(
        inputId = "stfr_select",  # ID of the 4th selectInput
        label = "Select a stfr column:",  # Label for the 4th selectInput
        choices = markers_stfr,  # Choices for the 4th selectInput
        selected = NULL  # Initial selected value for the 4th selectInput
      ),
      width = 3 #Width of the sidebar panel
    ),
    
    mainPanel(  # Contains the main panel outputs
      plotOutput(outputId = "olink_boxplot"),  # Output for the Olink boxplot plot
      textOutput(outputId = "olink_anova_result"),  # Output for the Olink ANOVA result
      plotOutput(outputId = "biochem_boxplot"),  # Output for the biochem boxplot plot
      textOutput(outputId = "biochem_anova_result"),  # Output for the biochem ANOVA result
      plotOutput(outputId = "hep_boxplot"),  # Output for the hepcidin boxplot plot
      textOutput(outputId = "hep_anova_result"),  # Output for the hepcidin ANOVA result
      plotOutput(outputId = "stfr_boxplot"),  # Output for the stfr boxplot plot
      textOutput(outputId = "stfr_anova_result")  # Output for the stfr ANOVA result
    )
  )
)
######################################### SERVER #############################################
server <- function(input, output) {
  # Function to select and filter data and run an ANOVA
  compute_model <- function(df, column_name) {
    
    filtered_df <- df %>%
      select(timepoint, ID, !!sym(column_name)) %>%
      filter(!is.na(.[[column_name]]) & is.finite(.[[column_name]]))  # Filter the data, removing rows with NA or infinite values
    
    # Perform linear mixed-effects regression and calculate ANOVA
    lmer_model <- lmer(formula = as.formula(paste(column_name, "~ timepoint + (1 | ID)")), data = filtered_df)
    anova_result <- anova(lmer_model, type = "III")
    p_value <- formatC(anova_result[["Pr(>F)"]][1], format = "f", digits = 5)  # Format the p-value
    p_value_display <- ifelse(as.numeric(p_value) < 0.05, paste0(p_value, "*"), p_value)  # Add asterisk if p-value is significant
    
    textual_result <- paste("Repeated Measures ANOVA p-value:", p_value_display)
    
    ## this put all the things in a list in case you might need it,
    # but if you only need the textual label, you can simply
    # return textual_result
    results <- list(
      lmer_model = lmer_model,
      anova_result = anova_result,
      p_value = p_value,
      p_value_display = p_value_display,
      textual_result = textual_result
    )
    
    return(results)
    
  }
  
  # Run the function on Olink data
  output$olink_anova_result <- renderText({
    selected_marker <- input$olink_select
    selected_column_name <- paste0(selected_marker, "_log10")
    model_stuff <- compute_model(olink_data, column_name = selected_column_name)
    return(model_stuff$textual_result)
    
    })
  
  # Run the function on Biochem data
  output$biochem_anova_result <- renderText({
    selected_marker <- input$biochem_select
    selected_column_name <- paste0(selected_marker, "_biochem")
    model_stuff <- compute_model(biochem_data, column_name = selected_column_name)
    return(model_stuff$textual_result)
    
  })  
  
  
  # Run the function on hep data
  output$hep_anova_result <- renderText({
    selected_marker <- input$hepcidin_select
    selected_column_name <- paste0(selected_marker, "_hep")
    model_stuff <- compute_model(hep_data, column_name = selected_column_name)
    return(model_stuff$textual_result)
    
  })  
  
  # Run the function on hep data
  output$stfr_anova_result <- renderText({
    selected_marker <- input$stfr_select
    selected_column_name <- paste0(selected_marker, "_stfr")
    model_stuff <- compute_model(stfr_data, column_name = selected_column_name)
    return(model_stuff$textual_result)
    
  })  
  
  # Function to make boxplot
  
  plot_boxplot <- function(df, column_name) {
    
    ggplot(data = df, aes(x = timepoint, y = .data[[column_name]], color = timepoint)) +
      geom_boxplot() +
      labs(x = "timepoint", y = column_name) +
      theme_classic()
  }
  
  # Run boxplot function on olink
  output$olink_boxplot <- renderPlot({
    selected_marker <- input$olink_select
    selected_column_name <- paste0(selected_marker, "_log10")
    plot_stuff <- plot_boxplot(olink_data, column_name = selected_column_name)
      print(plot_stuff) 
  })
  
  # Run boxplot function on biochem
  output$biochem_boxplot <- renderPlot({
    selected_marker <- input$biochem_select
    selected_column_name <- paste0(selected_marker, "_biochem")
    plot_stuff <- plot_boxplot(biochem_data, column_name = selected_column_name)
    print(plot_stuff) 
  })
  
  # Run boxplot function on hep
  output$hep_boxplot <- renderPlot({
    selected_marker <- input$hepcidin_select
    selected_column_name <- paste0(selected_marker, "_hep")
    plot_stuff <- plot_boxplot(hep_data, column_name = selected_column_name)
    print(plot_stuff) 
  })
  
  # Run boxplot function on stfr
  output$stfr_boxplot <- renderPlot({
    selected_marker <- input$stfr_select
    selected_column_name <- paste0(selected_marker, "_stfr")
    plot_stuff <- plot_boxplot(stfr_data, column_name = selected_column_name)
    print(plot_stuff) 
  })
  
}

  

### function 1###

  

######################################### Run the app! #############################################


shinyApp(ui = ui, server = server)



