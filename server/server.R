library(shiny)
library(ggplot2)
library(png)


shinyServer(function(input, output, session) {

    circuit <- eventReactive(input$execute, {
        circuit <- qc.parse_program(input$code)
        qc.compile(circuit)
        qc.run(circuit)

        return(circuit)
    })

    output$diagram <- renderImage({
        filename <- "circuit.png"

        qc <- circuit()
        qc.draw_to_file(qc, filename)

        image_size = dim(readPNG(filename))
        list(
            src = filename,
            height = min(session$clientData$output_diagram_height, image_size[1])
        )
    }, deleteFile = FALSE)

    output$prob_plot <- renderPlot({
        qc <- circuit()
        sv <- qc$statevector
    
        probs <- basis_probs(sv)
        state_names <- names(probs)
        state_probs <- unlist(probs, use.names=FALSE)
        
        dat <- data.frame(basis_states=state_names, probability=state_probs)

        ggplot(dat, aes(x=basis_states, y=probability)) + geom_col(fill=rgb(0.3, 0.6, 1.0)) + labs(x="Basis State", y="Probability") + ylim(0, 1) + coord_flip() 
    }, height = function() session$clientData$output_prob_plot_height)

    output$sv_plot <- renderPlot({
        qc <- circuit()
        sv <- qc$statevector

        amps <- basis_amps(sv)
        state_names <- names(amps)
        state_amps <- unlist(amps, use.names=FALSE)
        
        dat <- data.frame(basis_states=state_names, amplitude=state_amps)

        reals  <- Re(sv)
        imag   <- Im(sv)
        angles <- atan2(imag, reals)
        colors <- rad_to_rgb(angles)

        ggplot(dat, aes(x=basis_states, y=amplitude, fill=state_names)) + theme(legend.position='none') + geom_col() + labs(x="Basis State", y="Amplitude") + scale_fill_manual(values = colors) + ylim(0, 1)
    }, height = function() session$clientData$output_sv_plot_height)

    output$states <- renderPrint({
        qc <- circuit()
        round(qc$statevector, 3)
    })
})
