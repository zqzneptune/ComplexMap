/* ============================================================
   complexmap.js — Cytoscape.js Interactive Network
   ComplexMap v2.0
   ============================================================ */

(function () {
    "use strict";

    // ── Register extensions ──────────────────────────────────────────────────
    // Each block checks if the extension script is loaded before registering it.
    
    if (typeof cytoscape !== "undefined") {
        
        // 1. fCoSE (Fast Compound Spring Embedder)
        if (typeof cytoscapeFcose !== "undefined") {
            cytoscape.use(cytoscapeFcose);
        }

        // 2. Cola (Constraint-based Layout)
        // Note: 'cytoscapeCola' is the function from cytoscape-cola.js
        // It requires 'cola' (the engine) to be globally available.
        if (typeof cytoscapeCola !== "undefined") {
            cytoscape.use(cytoscapeCola);
        }

        // 3. SVG Export
        if (typeof cytoscapeSvg !== "undefined") {
            cytoscape.use(cytoscapeSvg);
        }
        
    } else {
        console.error("Cytoscape core not found. Extensions could not be registered.");
    }

    // ── Global state ─────────────────────────────────────────────────────────
    var cy;
    var currentLayout = "preset";
    var labelsVisible = true;

    // ── Helper: compute node min/max protein count for mapping ────────────────
    var elements = window.complexmapElements || [];
    var proteinCounts = elements
        .filter(function (el) { return el.data && !el.data.source; })
        .map(function (el) { return el.data.proteinCount || 1; });
    var minProt = Math.max(1, Math.min.apply(null, proteinCounts));
    var maxProt = Math.max(1, Math.max.apply(null, proteinCounts));

    function mapRange(val, inMin, inMax, outMin, outMax) {
        if (inMin === inMax) return (outMin + outMax) / 2;
        return ((val - inMin) / (inMax - inMin)) * (outMax - outMin) + outMin;
    }

    // ── Cytoscape Styles ───────────────────────────────────────────────────────
    var nodeStyle = {
        selector: "node",
        style: {
            "background-color": "data(colorHex)",
            "width": function (ele) {
                return mapRange(ele.data("proteinCount") || 1, minProt, maxProt, 18, 80);
            },
            "height": function (ele) {
                return mapRange(ele.data("proteinCount") || 1, minProt, maxProt, 18, 80);
            },
            "label": "data(label)",
            "color": "#ffffff",
            "font-size": 9,
            "text-valign": "center",
            "text-halign": "center",
            "text-wrap": "wrap",
            "text-max-width": "80px",
            "border-width": 1.5,
            "border-color": "rgba(255,255,255,0.25)",
            "transition-property": "opacity, background-color, width, height",
            "transition-duration": "200ms"
        }
    };

    var edgeStyle = {
        selector: "edge",
        style: {
            "width": function (ele) {
                return mapRange(ele.data("weight") || 0, 0, 1, 0.5, 5);
            },
            "line-color": "#556070",
            "opacity": function (ele) {
                return mapRange(ele.data("weight") || 0, 0, 1, 0.15, 0.9);
            },
            "curve-style": "bezier",
            "transition-property": "opacity, width",
            "transition-duration": "200ms"
        }
    };

    var selectedStyle = {
        selector: "node:selected",
        style: {
            "border-width": 3,
            "border-color": "#00d9ff",
            "background-color": "data(colorHex)"
        }
    };

    var fadedStyle = {
        selector: ".faded",
        style: { "opacity": 0.1 }
    };

    var highlightedStyle = {
        selector: ".highlighted",
        style: {
            "border-width": 3,
            "border-color": "#FFD700",
            "opacity": 1
        }
    };

    // ── Initialize Cytoscape ──────────────────────────────────────────────────
    document.addEventListener("DOMContentLoaded", function () {
        cy = cytoscape({
            container: document.getElementById("cy"),
            elements: elements,
            style: [nodeStyle, edgeStyle, selectedStyle, fadedStyle, highlightedStyle],
            layout: { name: "preset" },
            minZoom: 0.05,
            maxZoom: 8,
            wheelSensitivity: 0.3
        });

        // ── Node tap → send to Shiny ─────────────────────────────────────────
        cy.on("tap", "node", function (evt) {
            var node = evt.target;
            if (node.isParent && node.isParent()) return;
            if (typeof Shiny !== "undefined") {
                Shiny.setInputValue("selected_node", node.data("id"), { priority: "event" });
            }
        });

        // ── Click on canvas (deselect) ───────────────────────────────────────
        cy.on("tap", function (evt) {
            if (evt.target === cy) {
                clearHighlight();
                if (typeof Shiny !== "undefined") {
                    Shiny.setInputValue("selected_node", "", { priority: "event" });
                }
            }
        });

        cy.fit(undefined, 40);
    });

    // ── Layout Switcher ───────────────────────────────────────────────────────
    function applyLayout(name) {
        if (!cy) return;
        currentLayout = name;

        var layoutOptions = { 
            name: name, 
            animate: true, 
            animationDuration: 600 
        };

        if (name === "preset") {
            layoutOptions.positions = function (node) {
                return {
                    x: node.data("_preset_x") || node.position("x"),
                    y: node.data("_preset_y") || node.position("y")
                };
            };
        } 
        else if (name === "fcose") {
            layoutOptions.quality = "proof";
            layoutOptions.randomize = true;
            layoutOptions.nodeSeparation = 100;
        } 
        else if (name === "cola") {
            // Cola works best with specific spacing
            layoutOptions.nodeSpacing = 45;
            layoutOptions.edgeLength = 100;
            layoutOptions.maxSimulationTime = 2000;
        } 
        else if (name === "cose") {
            // Built-in CoSE requires higher repulsion for clear biological maps
            layoutOptions.nodeRepulsion = 400000;
            layoutOptions.idealEdgeLength = 100;
        }

        cy.layout(layoutOptions).run();
    }

    // ── Store preset positions after initial render ───────────────────────────
    document.addEventListener("DOMContentLoaded", function () {
        setTimeout(function () {
            if (!cy) return;
            cy.nodes().forEach(function (node) {
                node.data("_preset_x", node.position("x"));
                node.data("_preset_y", node.position("y"));
            });
        }, 300);
    });

    // ── Edge filtering ────────────────────────────────────────────────────────
    function applyEdgeFilter(minWeight) {
        if (!cy) return;
        cy.batch(function () {
            cy.edges().forEach(function (edge) {
                var w = parseFloat(edge.data("weight") || 0);
                if (w < minWeight) {
                    edge.style("display", "none");
                } else {
                    edge.style("display", "element");
                }
            });
        });
    }

    // ── Domain filter ─────────────────────────────────────────────────────────
    function applyDomainFilter(domain) {
        if (!cy) return;
        if (!domain || domain === "") {
            cy.nodes().removeClass("faded");
            cy.edges().removeClass("faded");
            return;
        }
        cy.nodes().each(function (node) {
            if (node.data("primaryFunctionalDomain") === domain) {
                node.removeClass("faded");
            } else {
                node.addClass("faded");
            }
        });
        cy.edges().each(function (edge) {
            if (edge.source().hasClass("faded") || edge.target().hasClass("faded")) {
                edge.addClass("faded");
            } else {
                edge.removeClass("faded");
            }
        });
    }

    // ── Protein search highlight ──────────────────────────────────────────────
    function highlightProtein(query) {
        if (!cy) return;
        clearHighlight();

        // Clear any previous search status
        if (typeof Shiny !== "undefined") {
            Shiny.setInputValue("protein_search_status", "", { priority: "event" });
        }

        if (!query || query.trim() === "") return;

        var q = query.trim().toLowerCase();

        var matchedNodes = cy.nodes().filter(function (node) {
            var proteins = (node.data("proteins") || "")
                .split(",")
                .map(function (p) { return p.trim().toLowerCase(); });

            var labelMatch = (node.data("label") || "").toLowerCase() === q;

            return proteins.indexOf(q) !== -1 || labelMatch;
        });

        if (matchedNodes.length === 0) {
            if (typeof Shiny !== "undefined") {
                Shiny.setInputValue("protein_search_status", "No hits found for: " + query.trim(), { priority: "event" });
            }
            return;
        }

        // Report how many nodes matched
        if (typeof Shiny !== "undefined") {
            Shiny.setInputValue("protein_search_status", 
                matchedNodes.length + " complex" + (matchedNodes.length > 1 ? "es" : "") + " found",
                { priority: "event" });
        }

        cy.nodes().addClass("faded");
        cy.edges().addClass("faded");
        matchedNodes.removeClass("faded").addClass("highlighted");
        matchedNodes.connectedEdges().removeClass("faded");
    }

    function clearHighlight() {
        if (!cy) return;
        cy.nodes().removeClass("faded highlighted");
        cy.edges().removeClass("faded");
    }

    // ── Label toggle ──────────────────────────────────────────────────────────
    // FIX: removed syncLabelState() entirely — initial state is now pushed by
    // R via session$onFlushed, which fires only after Cytoscape is ready.
    // toggleLabels() is the single source of truth for label visibility.
    function toggleLabels(show) {
        if (!cy) return;
        labelsVisible = show;
        cy.nodes().style("label", show ? "data(label)" : "");
        console.log("Label toggle called:", show ? "showing" : "hiding", "labels");
    }

    // ── PNG / SVG export ──────────────────────────────────────────────────────
    function exportPNG() {
        if (!cy) {
            console.error("Cytoscape instance not available");
            if (typeof Shiny !== "undefined") {
                Shiny.setInputValue("export_status", "PNG export failed: Cytoscape not ready", { priority: "event" });
            }
            return;
        }
        try {
            var dataUrl = cy.png({ full: true, scale: 2, bg: "#1a1d2e" });
            var link = document.createElement("a");
            link.download = "complexmap.png";
            link.href = dataUrl;
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            if (typeof Shiny !== "undefined") {
                Shiny.setInputValue("export_status", "PNG exported successfully", { priority: "event" });
            }
        } catch (e) {
            console.error("Export PNG failed:", e);
            if (typeof Shiny !== "undefined") {
                Shiny.setInputValue("export_status", "PNG export failed: " + e.message, { priority: "event" });
            }
        }
    }

    function exportSVG() {
        if (!cy) {
            console.error("Cytoscape instance not available");
            if (typeof Shiny !== "undefined") {
                Shiny.setInputValue("export_status", "SVG export failed: Cytoscape not ready", { priority: "event" });
            }
            return;
        }

        // Re-attempt extension registration in case it wasn't ready at init time
        if (typeof cy.svg !== "function" && typeof cytoscapeSvg !== "undefined") {
            try { cytoscape.use(cytoscapeSvg); } catch (e) { /* already registered */ }
        }

        if (typeof cy.svg !== "function") {
            console.error("SVG export not available. Check if cytoscape-svg is loaded.");
            if (typeof Shiny !== "undefined") {
                Shiny.setInputValue("export_status", "SVG export failed: cytoscape-svg extension not loaded", { priority: "event" });
            }
            return;
        }

        try {
            var svgContent = cy.svg({ full: true, bg: "#1a1d2e" });
            var blob = new Blob([svgContent], { type: "image/svg+xml;charset=utf-8" });
            var url = URL.createObjectURL(blob);
            var link = document.createElement("a");
            link.download = "complexmap.svg";
            link.href = url;
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            setTimeout(function () { URL.revokeObjectURL(url); }, 100);
            if (typeof Shiny !== "undefined") {
                Shiny.setInputValue("export_status", "SVG exported successfully", { priority: "event" });
            }
        } catch (e) {
            console.error("Export SVG failed:", e);
            if (typeof Shiny !== "undefined") {
                Shiny.setInputValue("export_status", "SVG export failed: " + e.message, { priority: "event" });
            }
        }
    }

    // ── Shiny Message Handlers ────────────────────────────────────────────────
    if (typeof Shiny !== "undefined") {
        Shiny.addCustomMessageHandler("change_layout", function (msg) {
            applyLayout(msg.name);
        });

        Shiny.addCustomMessageHandler("filter_edges", function (msg) {
            applyEdgeFilter(msg.minWeight || 0);
        });

        Shiny.addCustomMessageHandler("filter_domain", function (msg) {
            applyDomainFilter(msg.domain);
        });

        Shiny.addCustomMessageHandler("highlight_protein", function (msg) {
            highlightProtein(msg.query);
        });

        // Single handler for both initial state and subsequent toggles
        Shiny.addCustomMessageHandler("toggle_labels", function (msg) {
            toggleLabels(msg.show);
        });

        Shiny.addCustomMessageHandler("export_png", function (msg) {
            exportPNG();
        });

        Shiny.addCustomMessageHandler("export_svg", function (msg) {
            exportSVG();
        });
    }

})();
