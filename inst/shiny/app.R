## inst/shiny/app.R — ComplexMap Interactive Explorer (Cytoscape.js)
library(shiny)
library(jsonlite)
library(htmltools)

# Retrieve data stored by runComplexMapApp() via shinyOptions
elements_json  <- getShinyOption("complexmap_elements",  default = "[]")
nodes_df       <- getShinyOption("complexmap_nodes",     default = data.frame())
edges_df       <- getShinyOption("complexmap_edges",     default = data.frame())
domains        <- getShinyOption("complexmap_domains",   default = character(0))
layout_info    <- getShinyOption("complexmap_layout_info", default = list())

# Compute slider ranges from edge data
max_weight     <- if (nrow(edges_df) > 0) max(edges_df$weight, na.rm = TRUE) else 1
min_weight     <- if (nrow(edges_df) > 0) min(edges_df$weight, na.rm = TRUE) else 0

# ──────────────────────────────────────────────────────────────────────────────
# UI
# ──────────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  title = "ComplexMap Explorer",

  # CSS & JS dependencies
  tags$head(
    tags$link(rel = "stylesheet", href = "style.css"),
    tags$script(src = "cytoscape.min.js"),
    tags$script(src = "layout-base.js"),
    tags$script(src = "cose-base.js"),
    tags$script(src = "cytoscape-fcose.min.js"),
    tags$script(src = "cytoscape-svg.min.js"),

    # Inject network data before complexmap.js loads
    tags$script(
      HTML(paste0("window.complexmapElements = ", elements_json, ";"))
    ),
    # Inject unique protein names for autocomplete
    tags$script(HTML(paste0(
      "window.complexmapProteins = ",
      jsonlite::toJSON(sort(unique(unlist(strsplit(nodes_df$proteins, ",")))), auto_unbox = FALSE),
      ";"
    ))),
    tags$script(src = "complexmap.js")
  ),

  # ── Layout: sidebar + main
  div(class = "cm-wrapper",

    # ── Sidebar ──────────────────────────────────────────────────────────────
    div(class = "cm-sidebar",
      div(class = "cm-logo", "ComplexMap", tags$span("Explorer")),

      # Layout selector
      div(class = "cm-section",
        tags$label("Layout", class = "cm-label"),
        selectInput("layout_choice", label = NULL,
          choices  = c("Preset (FR)" = "preset", "fCoSE" = "fcose"),
          selected = "preset"
        )
      ),

      # Edge weight filter
      div(class = "cm-section",
        tags$label("Min. Edge Weight", class = "cm-label"),
        sliderInput("edge_weight_min", label = NULL,
          min = floor(min_weight * 10) / 10,
          max = ceiling(max_weight * 10) / 10,
          value = floor(min_weight * 10) / 10,
          step = 0.01
        )
      ),

      # Functional domain filter
      if (length(domains) > 0) {
        div(class = "cm-section",
          tags$label("Functional Domain", class = "cm-label"),
          selectInput("domain_filter", label = NULL,
            choices = c("All" = "", domains), selected = ""
          )
        )
      },

      # Protein search
      div(class = "cm-section",
          tags$label("Search Protein", class = "cm-label"),
          div(class = "cm-search-wrapper",
              div(class = "cm-search-row",
                  div(style = "position: relative; flex: 1;",
                      textInput("protein_search", label = NULL, placeholder = "e.g. TP53"),
                      tags$ul(id = "protein_suggestions", class = "cm-suggestions")
                  )
              ),
              uiOutput("protein_search_status")
          ),
          tags$script(HTML("
        (function() {
          var allProteins = window.complexmapProteins || [];
          var input      = null;
          var list       = null;
          var activeIdx  = -1;
    
          function init() {
            input = document.getElementById('protein_search');
            list  = document.getElementById('protein_suggestions');
            if (!input || !list) return;
    
            // Show suggestions on typing
            input.addEventListener('input', function() {
              var q = input.value.trim().toLowerCase();
              renderSuggestions(q);
              activeIdx = -1;
            });
    
            // Keyboard navigation: ArrowDown / ArrowUp / Enter / Escape
            input.addEventListener('keydown', function(e) {
              var items = list.querySelectorAll('li');
              if (e.key === 'ArrowDown') {
                e.preventDefault();
                activeIdx = Math.min(activeIdx + 1, items.length - 1);
                updateActive(items);
              } else if (e.key === 'ArrowUp') {
                e.preventDefault();
                activeIdx = Math.max(activeIdx - 1, 0);
                updateActive(items);
              } else if (e.key === 'Enter') {
                if (activeIdx >= 0 && items[activeIdx]) {
                  // Select the highlighted suggestion
                  selectSuggestion(items[activeIdx].textContent);
                } else {
                  // No suggestion selected — run search directly
                  hideSuggestions();
                  Shiny.setInputValue('protein_search_enter', Date.now(), { priority: 'event' });
                }
              } else if (e.key === 'Escape') {
                hideSuggestions();
              }
            });
    
            // Hide suggestions when clicking outside
            document.addEventListener('click', function(e) {
              if (!input.contains(e.target) && !list.contains(e.target)) {
                hideSuggestions();
              }
            });
          }
    
          function renderSuggestions(q) {
            list.innerHTML = '';
            if (!q) { hideSuggestions(); return; }
    
            // Match from beginning of gene name, case-insensitive
            var matches = allProteins.filter(function(p) {
              return p.toLowerCase().indexOf(q) === 0;
            }).slice(0, 10); // cap at 10 suggestions
    
            if (matches.length === 0) { hideSuggestions(); return; }
    
            matches.forEach(function(protein) {
              var li = document.createElement('li');
              // Bold the matched prefix, preserve original casing
              li.innerHTML = '<strong>' + protein.slice(0, q.length) + '</strong>' + protein.slice(q.length);
              li.addEventListener('mousedown', function(e) {
                e.preventDefault(); // prevent blur before click registers
                selectSuggestion(protein);
              });
              list.appendChild(li);
            });
    
            list.style.display = 'block';
          }
    
          function selectSuggestion(protein) {
            // Push value into Shiny input and trigger search
            input.value = protein;
            Shiny.setInputValue('protein_search', protein, { priority: 'event' });
            hideSuggestions();
            Shiny.setInputValue('protein_search_enter', Date.now(), { priority: 'event' });
          }
    
          function updateActive(items) {
            items.forEach(function(li, i) {
              li.classList.toggle('active', i === activeIdx);
            });
          }
    
          function hideSuggestions() {
            if (list) list.style.display = 'none';
            activeIdx = -1;
          }
    
          // Init after DOM is ready
          if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', init);
          } else {
            init();
          }
        })();
      "))
      ),

      # Label toggle
      div(class = "cm-section",
        checkboxInput("show_labels", "Show node labels", value = TRUE)
      ),

      tags$hr(class = "cm-divider"),

      # Export buttons
      div(class = "cm-section",
        tags$label("Export", class = "cm-label"),
        div(class = "cm-btn-group",
          actionButton("export_png", "PNG", class = "cm-btn cm-btn-export"),
          actionButton("export_svg", "SVG", class = "cm-btn cm-btn-export"),
          downloadButton("export_json", "JSON", class = "cm-btn cm-btn-export"),
          downloadButton("export_tsv",  "TSV",  class = "cm-btn cm-btn-export")
        )
      ),

      # Export status message
      div(class = "cm-section",
        uiOutput("export_status")
      ),

      # Layout info
      if (length(layout_info) > 0) {
        div(class = "cm-info-footer",
          sprintf("Layout: %s | Seed: %s",
                  layout_info$method %||% "fr",
                  layout_info$seed %||% "123")
        )
      }
    ), # end sidebar

    # ── Main canvas ──────────────────────────────────────────────────────────
    div(class = "cm-main",
      div(id = "cy"),

      # Node detail panel (appears on node click)
      conditionalPanel(
        condition = "input.selected_node != ''",
        div(class = "cm-detail-panel",
          uiOutput("node_detail")
        )
      )
    )
  ) # end wrapper
)

# ──────────────────────────────────────────────────────────────────────────────
# Null-coalescing helper
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ──────────────────────────────────────────────────────────────────────────────
# Server
# ──────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Send layout command to JS ───────────────────────────────────────────────
  observeEvent(input$layout_choice, {
    session$sendCustomMessage("change_layout", list(name = input$layout_choice))
  })

  # ── Edge weight filter → JS ─────────────────────────────────────────────────
  observeEvent(input$edge_weight_min, {
    session$sendCustomMessage("filter_edges", list(minWeight = input$edge_weight_min))
  })

  # ── Domain filter → JS ─────────────────────────────────────────────────────
  observeEvent(input$domain_filter, {
    session$sendCustomMessage("filter_domain", list(domain = input$domain_filter))
  })

  # ── Protein search → JS ────────────────────────────────────────────────────
  # FIX: also trigger on Enter key in the text box (input$protein_search change),
  # not only on the Go button, and pass trimmed query.
  # Shared search logic
  doSearch <- function() {
    query <- trimws(input$protein_search)
    session$sendCustomMessage("highlight_protein", list(query = query))
  }


  # Trigger on Enter key (fires when textInput value changes and focus is committed)
  # Using a dedicated JS binding is the most reliable approach for Enter key detection
  observeEvent(input$protein_search_enter, {
    doSearch()
  })

  # Clear highlights when search box is emptied
  observeEvent(input$protein_search, {
    if (trimws(input$protein_search) == "") {
      session$sendCustomMessage("highlight_protein", list(query = ""))
    }
  }, ignoreInit = TRUE)

  output$protein_search_status <- renderUI({
      status <- input$protein_search_status
      if (is.null(status) || status == "") return(NULL)

      is_error <- grepl("No hits", status)
      div(
          class = if (is_error) "cm-status cm-status-error" else "cm-status cm-status-success",
          status
      )
  })
  # ── Label toggle → JS ──────────────────────────────────────────────────────
  # FIX: Use session$onFlushed to ensure the first message is sent only after
  # Shiny has finished its initial render and the JS/Cytoscape instance is ready.
  # Subsequent toggles fire immediately via observeEvent as normal.
  
  # session$onFlushed must be called at the top level of server, NOT inside
  # observe() or any other reactive context. once = TRUE fires only after the
  # first complete render, guaranteeing Cytoscape is mounted before we send.
  session$onFlushed(function() {
    session$sendCustomMessage("toggle_labels", list(show = TRUE))
  }, once = TRUE)

  # All subsequent user-driven toggles — ignoreInit skips the startup firing
  observeEvent(input$show_labels, {
    session$sendCustomMessage("toggle_labels", list(show = isTRUE(input$show_labels)))
  }, ignoreInit = TRUE)

  # ── Node click: display detail panel ───────────────────────────────────────
  output$node_detail <- renderUI({
    node_id <- input$selected_node
    if (is.null(node_id) || node_id == "") return(NULL)

    nd <- nodes_df[nodes_df$complexId == node_id, ]
    if (nrow(nd) == 0) return(div("Node not found."))

    proteins_list <- strsplit(nd$proteins[1], ",")[[1]]

    tagList(
      tags$button(class = "cm-close-btn",
                  onclick = "Shiny.setInputValue('selected_node', '')",
                  HTML("&times;")),
      div(class = "cm-detail-title", node_id),
      div(class = "cm-detail-row",
        span(class = "cm-detail-key", "Function:"),
        span(class = "cm-detail-val", nd$primaryFunctionalDomain[1])
      ),
      div(class = "cm-detail-row",
        span(class = "cm-detail-key", "Protein Count:"),
        span(class = "cm-detail-val", nd$proteinCount[1])
      ),
      if ("betweenness" %in% names(nd)) {
        div(class = "cm-detail-row",
          span(class = "cm-detail-key", "Betweenness:"),
          span(class = "cm-detail-val", round(nd$betweenness[1], 4))
        )
      },
      if ("degree" %in% names(nd)) {
        div(class = "cm-detail-row",
          span(class = "cm-detail-key", "Degree:"),
          span(class = "cm-detail-val", nd$degree[1])
        )
      },
      if ("topEnrichedFunctions" %in% names(nd) && !is.na(nd$topEnrichedFunctions[1])) {
        tagList(
          div(class = "cm-detail-subheader", "Top Enriched Functions"),
          div(class = "cm-detail-enrichment",
              nd$topEnrichedFunctions[1])
        )
      },
      div(class = "cm-detail-subheader",
          paste0("Members (", length(proteins_list), ")")),
      div(class = "cm-proteins-list",
          tags$ul(lapply(proteins_list, tags$li))
      )
    )
  })

  # ── Export PNG / SVG (triggered by JS directly) ────────────────────────────
  observeEvent(input$export_png, {
    session$sendCustomMessage("export_png", list())
  })

  observeEvent(input$export_svg, {
    session$sendCustomMessage("export_svg", list())
  })

  # ── Export status display ──────────────────────────────────────────────────
  output$export_status <- renderUI({
    if (is.null(input$export_status) || input$export_status == "") return(NULL)

    status_class <- if (grepl("failed", input$export_status, ignore.case = TRUE))
                    "cm-status-error" else "cm-status-success"

    div(class = paste0("cm-status ", status_class),
        input$export_status)
  })

  # ── Export JSON ────────────────────────────────────────────────────────────
  output$export_json <- downloadHandler(
    filename = function() paste0("complexmap_", Sys.Date(), ".json"),
    content  = function(file) writeLines(elements_json, file)
  )

  # ── Export TSV ─────────────────────────────────────────────────────────────
  output$export_tsv <- downloadHandler(
    filename = function() paste0("complexmap_nodes_", Sys.Date(), ".tsv"),
    content  = function(file) {
      utils::write.table(nodes_df, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
}

shinyApp(ui, server)
