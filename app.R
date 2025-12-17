# Required Libraries
library(shiny)
library(rentrez)
library(xml2)
library(tibble)
library(dplyr)
library(openxlsx)
library(shinythemes)

# Helper function: Build PubMed query
build_query <- function(search_term, additional_terms = NULL) {
  # Start with main term in MeSH
  query_parts <- paste0(search_term, "[MeSH Terms]")
  
  # Add main term in Title/Abstract
  query_parts <- c(query_parts, paste0(search_term, "[Title/Abstract]"))
  
  # Add additional terms if provided
  if (!is.null(additional_terms) && length(additional_terms) > 0) {
    additional_queries <- paste0(additional_terms, "[Title/Abstract]")
    query_parts <- c(query_parts, additional_queries)
  }
  
  # Combine with OR
  paste0("(", paste(query_parts, collapse = " OR "), ")")
}

# Function 1: Fetch and parse to tibble
fetch_pubmed_data <- function(search_term, additional_terms = NULL, retmax = 300) {
  tryCatch({
    query <- build_query(search_term, additional_terms)
    
    # Cap at 300 to prevent crashes
    retmax <- min(retmax, 300)
    
    search <- entrez_search(
      db = "pubmed",
      term = query,
      retmax = retmax
    )
    
    if (length(search$ids) == 0) {
      warning("No records found for ", search_term)
      return(NULL)
    }
    
    # Respect rate limits (3 requests/second without API key)
    Sys.sleep(0.34)
    
    xml_raw <- entrez_fetch(
      db = "pubmed",
      id = search$ids,
      rettype = "abstract",
      retmode = "xml"
    )
    
    xml_doc <- read_xml(xml_raw)
    articles <- xml_find_all(xml_doc, ".//PubmedArticle")
    
    result <- tibble(
      Search_Term = search_term,
      PMID = xml_text(xml_find_first(articles, ".//PMID")),
      Title = xml_text(xml_find_first(articles, ".//ArticleTitle")),
      Abstract = sapply(articles, function(a) {
        abstract_parts <- xml_find_all(a, ".//AbstractText")
        if (length(abstract_parts) == 0) return(NA_character_)
        paste(xml_text(abstract_parts), collapse = " ")
      }),
      Journal = xml_text(xml_find_first(articles, ".//Journal/Title")),
      Year = xml_text(xml_find_first(articles, ".//PubDate/Year")),
      DOI = sapply(articles, function(a) {
        doi_node <- xml_find_first(a, ".//ArticleId[@IdType='doi']")
        if (is.na(doi_node)) return(NA_character_)
        xml_text(doi_node)
      }),
      Authors = sapply(articles, function(a) {
        authors <- xml_find_all(a, ".//Author")
        if (length(authors) == 0) return(NA_character_)
        author_names <- sapply(authors, function(auth) {
          last <- xml_text(xml_find_first(auth, ".//LastName"))
          first <- xml_text(xml_find_first(auth, ".//ForeName"))
          if (is.na(last)) return(NA_character_)
          if (is.na(first)) return(last)
          paste(last, first, sep = ", ")
        })
        paste(author_names[!is.na(author_names)], collapse = "; ")
      })
    )
    
    return(result)
    
  }, error = function(e) {
    warning("Error fetching data for ", search_term, ": ", e$message)
    return(NULL)
  })
}

# Function 2: Save raw XML file (to temp for download)
save_pubmed_xml <- function(search_term, additional_terms = NULL, retmax = 300) {
  tryCatch({
    query <- build_query(search_term, additional_terms)
    
    # Cap at 300 to prevent crashes
    retmax <- min(retmax, 300)
    
    search <- entrez_search(
      db = "pubmed",
      term = query,
      retmax = retmax
    )
    
    if (length(search$ids) == 0) {
      warning("No records found for ", search_term)
      return(NULL)
    }
    
    Sys.sleep(0.34)
    
    xml_raw <- entrez_fetch(
      db = "pubmed",
      id = search$ids,
      rettype = "abstract",
      retmode = "xml"
    )
    
    # Create PostgreSQL-friendly XML wrapper
    xml_wrapped <- paste0(
      '<?xml version="1.0" encoding="UTF-8"?>\n',
      '<pubmed_export>\n',
      '  <metadata>\n',
      '    <search_term>', search_term, '</search_term>\n',
      '    <export_date>', Sys.Date(), '</export_date>\n',
      '    <record_count>', length(search$ids), '</record_count>\n',
      '  </metadata>\n',
      '  <records>\n',
      xml_raw,
      '\n  </records>\n',
      '</pubmed_export>'
    )
    
    # Save to temp file for download
    outfile <- tempfile(pattern = paste0(gsub("[^A-Za-z0-9]", "_", search_term), "_pubmed_"), fileext = ".xml")
    writeLines(xml_wrapped, outfile, useBytes = TRUE)
    
    message("Created XML for ", search_term, " (", length(search$ids), " records)")
    return(outfile)
    
  }, error = function(e) {
    warning("Error saving XML for ", search_term, ": ", e$message)
    return(NULL)
  })
}

# Shiny UI
ui <- fluidPage(
  theme = shinytheme("lumen"),
  titlePanel("PubMed Crawler 🕵️"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("search_term", "Search Term (required):", value = "ZEB2", 
                placeholder = "e.g., ZEB2, Parkinson's disease, CRISPR"),
      textInput("additional_terms", "Additional Terms (optional, comma-separated):", 
                value = "",
                placeholder = "e.g., SIP1, synonyms, related terms"),
      numericInput("retmax", "Maximum Results:", value = 300, 
                   min = 1, max = 300, step = 50),
      helpText("Note: Maximum 300 papers to prevent app crashes."),
      hr(),
      actionButton("fetch", "Fetch Papers", class = "btn-primary"),
      br(), br(),
      downloadButton("download_xlsx", "Download XLSX"),
      downloadButton("download_xml", "Download XML"),
      hr(),
      helpText("Search for genes, diseases, cell types, methods,..."),
      helpText("Large queries may take time."),
      hr(),
      h4("Andrea Conidi - 2025"),
      hr(),
      h4("github:")
    ),
    
    mainPanel(
      h4("Results Summary"),
      verbatimTextOutput("summary"),
      hr(),
      h4("Preview (First 10 Results)"),
      tableOutput("preview")
     
    )
  )
)

# Shiny Server
server <- function(input, output, session) {
  
  # Reactive values to store results
  results <- reactiveValues(
    data = NULL,
    xml_file = NULL
  )
  
  # Fetch data when button clicked
  observeEvent(input$fetch, {
    req(input$search_term)  # Search term is required
    
    # Validate search term input
    if (trimws(input$search_term) == "") {
      showModal(modalDialog(
        title = "Error",
        "Search term is required!",
        easyClose = TRUE
      ))
      return()
    }
    
    withProgress(message = 'Fetching from PubMed...', value = 0, {
      
      # Parse additional terms (may be empty)
      additional <- NULL
      if (!is.null(input$additional_terms) && trimws(input$additional_terms) != "") {
        additional <- trimws(unlist(strsplit(input$additional_terms, ",")))
        additional <- additional[additional != ""]  # Remove empty strings
        if (length(additional) == 0) additional <- NULL
      }
      
      # Fetch data
      incProgress(0.3, detail = "Retrieving papers...")
      data <- fetch_pubmed_data(input$search_term, additional, input$retmax)
      
      # Save XML
      incProgress(0.6, detail = "Generating XML...")
      xml_file <- save_pubmed_xml(input$search_term, additional, input$retmax)
      
      # Store results
      results$data <- data
      results$xml_file <- xml_file
      
      incProgress(1, detail = "Complete!")
    })
  })
  
  # Summary output
  output$summary <- renderText({
    if (is.null(results$data)) {
      return("No data fetched yet. Enter a search term and click 'Fetch Papers' to start.")
    }
    
    total_found <- nrow(results$data)
    capped_msg <- if (input$retmax >= 300) {
      "\n(Limited to 300 papers maximum)"
    } else {
      ""
    }
    
    paste0(
      "Search Term: ", input$search_term, "\n",
      "Total Papers Retrieved: ", total_found, capped_msg, "\n",
      "Papers with Abstracts: ", sum(!is.na(results$data$Abstract)), "\n",
      "Date Range: ", min(results$data$Year, na.rm = TRUE), " - ",
      max(results$data$Year, na.rm = TRUE)
    )
  })
  
  # Preview table
  output$preview <- renderTable({
    req(results$data)
    head(results$data %>% select(PMID, Title, Year, Journal), 10)
  })
  
  # Download XLSX
  output$download_xlsx <- downloadHandler(
    filename = function() {
      clean_term <- gsub("[^A-Za-z0-9]", "_", input$search_term)
      paste0(clean_term, "_pubmed_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      req(results$data)
      write.xlsx(results$data, file)
    }
  )
  
  # Download XML
  output$download_xml <- downloadHandler(
    filename = function() {
      clean_term <- gsub("[^A-Za-z0-9]", "_", input$search_term)
      paste0(clean_term, "_pubmed_", Sys.Date(), ".xml")
    },
    content = function(file) {
      req(results$xml_file)
      file.copy(results$xml_file, file)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)