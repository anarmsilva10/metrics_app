from shiny import App, ui, module, reactive, render, req
from shinywidgets import output_widget, render_widget
import shinyswatch
import pandas as pd
import io
import plotly.graph_objects as go
import atexit
from backend.depth import filter_bed
from backend.coverage_metrics import calculate_metrics
from backend.s3 import list_samples, list_files, download_files
from backend.single_gene import get_genes, get_exons, create_bed
from backend.gene_panel import get_gene_panel, get_panel_genes, get_panel_exons, create_panel_bed
from backend.plot import coverage_plot
from backend.helper import clean_temp_dirs, update_file
from log import shared_log

# Global variables
TITLE = "NGS Metrics Calculator"
NS = "metrics_app"

# Logging initiation
shared_log.logger = shared_log.Logging(basename='metrics_app', foldername='logs', log_level='INFO').get_logger()

update_file()

@module.ui
def metrics_app_ui(show_title=False):
    icon = ui.HTML('''<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-chevron-double-down" viewBox="0 0 16 16">
    <path fill-rule="evenodd" d="M1.646 6.646a.5.5 0 0 1 .708 0L8 12.293l5.646-5.647a.5.5 0 0 1 .708.708l-6 6a.5.5 0 0 1-.708 0l-6-6a.5.5 0 0 1 0-.708"/>
    <path fill-rule="evenodd" d="M1.646 2.646a.5.5 0 0 1 .708 0L8 8.293l5.646-5.647a.5.5 0 0 1 .708.708l-6 6a.5.5 0 0 1-.708 0l-6-6a.5.5 0 0 1 0-.708"/>
    </svg>''')

    popover_style = ui.tags.style(""".metrics {
                                max-height: 200px;
                                overflow-y: auto;
                                padding-right:10px}""")

    # Somatics Cases Page
    page_somatics = ui.layout_column_wrap(
                # Input Card
                ui.card(
                    ui.card_header("Input Panel:"),
                        # Input Files
                        ui.input_file('cram_file', 'Upload BAM/CRAM File'),
                        ui.input_file('bed_file', 'Upload BED File'),
                        # Select Method
                        ui.input_select(
                            'method',
                            'Choose Depth Extraction Method',
                            {
                                'run_pysam_coverage': 'Pysam Coverage',
                                'run_pysam_depth': 'Pysam depth',
                                'run_samtools_depth': 'Samtools depth'
                            }
                        ),
                        ui.input_switch("stats_somatics", "Statistics Metrics", True),
                        # Button to run backend scripts to calculate metrics
                        ui.input_action_button('run', label='Run Metrics Analysis'),
                        # Download output
                        ui.input_text('output_somatics', 'Output file', placeholder='e.g. metrics.xlsx'),
                        ui.download_button('download_data_somatics', 'Download Report'),
                    height='auto'
                ),
                # Metrics Display Card
                ui.card(
                    ui.card_header('Results Panel:'),
                    popover_style,
                    ui.popover(
                        ui.span(icon, 'Metrics Selection', style='position:absolute; top: 10px; right: 7px; font-size:12px;'),
                        ui.output_ui('metrics_select_somatics'),
                        placement='right'),
                    ui.output_data_frame('display_output_somatics'),
                    height='auto'
                )
            )

    # WES Cases Page
    page_WES = ui.navset_card_underline(
        ui.nav_panel("Single Gene",
                    ui.layout_column_wrap(
                        # Input Card
                        ui.card(
                            ui.card_header("Input Panel:"),
                                # Input Files
                                ui.output_ui('update_samples'),
                                ui.output_ui('update_files'),
                                ui.input_radio_buttons('genome', 'Choose Reference Genome:',
                                                    {
                                                        'hg19': 'Genome Assembly GRCh37',
                                                        'hg38': 'Genome Assembly GRCh38'
                                                    },
                                                    # In order to not have an already selected genome and both options aligned
                                                    selected="", inline=True),
                                ui.input_selectize('gene_select', 'Choose Gene', choices=[]),
                                ui.input_switch('trim_utr_wes', "Trimming 5' and 3' UTR", False),
                                ui.input_slider('padding_wes', "Choose Padding:", 0, 20, 10),
                                ui.output_ui('exon_selector'),
                                ui.input_checkbox('select_exons', 'Select/Deselect all exons', False),
                                #Select Method
                                ui.input_select(
                                    'method_wes',
                                    'Choose Depth Extraction Method',
                                    {
                                        'run_pysam_coverage': 'Pysam Coverage',
                                        'run_pysam_depth': 'Pysam depth',
                                        'run_samtools_depth': 'Samtools depth'
                                    }
                                ),
                                ui.input_switch("stats_wes", "Statistics Metrics", False),
                                # Button to run backend scripts to calculate metrics
                                ui.input_action_button('run_wes', label='Run Metrics Analysis'),
                                # Plot
                                ui.input_checkbox('plot_wes', 'Display Coverage Plot', False),
                                # Download output
                                ui.input_text('output_wes', 'Output file', placeholder='e.g. metrics.xlsx'),
                                ui.download_button('download_data_wes', 'Download Report'),
                            height='auto'
                        ),
                        # Metrics Display Card
                        ui.card(
                            ui.card_header('Results Panel:'),
                            ui.card(
                                ui.card_header('Metrics Display:'),
                                popover_style,
                                ui.popover(
                                    ui.span(icon, 'Metrics Selection', style='position:absolute; top: 10px; right: 7px; font-size:12px;'),
                                    ui.output_ui('metrics_select_wes'),
                                    placement='right'
                                ),
                                ui.output_data_frame('display_output_wes'),
                            ),
                            ui.panel_conditional(
                                "input.plot_wes",
                                ui.card(
                                    ui.card_header('Coverage Plot:'),
                                    ui.input_slider("threshold_wes", 'Choose Coverage Threshold:', 1, 400, 20),
                                    ui.input_selectize('plot_exons_wes', 'Select Exons to Plot:', choices=[], selected=[], multiple=True),
                                    output_widget("plot_coverage_wes", width="100%", height="700px")
                                )
                            ), height='auto'
                    ))
                ),
        ui.nav_panel("Gene Panel",
                    ui.layout_column_wrap(
                        # Input Card
                        ui.card(
                            ui.card_header("Input Panel:"),
                                # Input Files
                                ui.output_ui('update_samples_gp'),
                                ui.output_ui('update_files_gp'),
                                ui.input_radio_buttons('genome_gp', 'Choose Reference Genome:',
                                                    {
                                                        'hg19': 'Genome Assembly GRCh37',
                                                        'hg38': 'Genome Assembly GRCh38'
                                                    },
                                                    # In order to not have an already selected genome and both options aligned
                                                    selected="", inline=True),
                                ui.output_ui('gene_panel'),
                                ui.input_switch('trim_utr_gp', "Trimming 5' and 3' UTR", False),
                                ui.input_slider('padding_gp', "Choose Padding:", 0, 20, 10),
                                # Select Method
                                ui.input_select(
                                    'method_gp',
                                    'Choose Depth Extraction Method',
                                    {
                                        'run_pysam_coverage': 'Pysam Coverage',
                                        'run_pysam_depth': 'Pysam depth',
                                        'run_samtools_depth': 'Samtools depth'
                                    }
                                ),
                                ui.input_switch("stats_gp", "Statistics Metrics", False),
                                # Gene and Exon Selection
                                ui.output_ui('select_gene'),
                                ui.output_ui('exon_select'),
                                ui.input_checkbox('select_exons_gp', 'Select/Deselect all exons', False),
                                # Button to run backend scripts to calculate metrics
                                ui.input_action_button('run_gp', label='Run Metrics Analysis'),
                                # Plot
                                ui.input_checkbox('plot_gp', 'Display Coverage Plot', False),
                                # Download output
                                ui.input_text('output_gp', 'Output file', placeholder='e.g. metrics.xlsx'),
                                ui.download_button('download_data_gp', 'Download Report'),
                            height='auto'
                        ),
                        # Metrics Display Card
                        ui.card(
                            ui.card_header('Results Panel:'),
                            ui.card(
                                ui.card_header('Metrics Display:'),
                                popover_style,
                                ui.popover(
                                    ui.span(icon, 'Metrics Selection', style='position:absolute; top: 10px; right: 7px; font-size:12px;'),
                                    ui.output_ui('metrics_select_gp')
                                ),
                                ui.output_data_frame('display_output_gp')
                            ),
                            ui.panel_conditional(
                                "input.plot_gp",
                                ui.card(
                                    ui.card_header('Coverage Plot:'),
                                    ui.input_slider("threshold_gp", 'Choose Coverage Threshold:', 1, 400, 20),
                                    ui.output_ui('gene_select_plotly_gp'),
                                    ui.output_ui('plot_exons_selector_gp'),
                                    output_widget("plot_coverage_gp", width="100%", height="100%")
                                )
                            ), height='auto'
                    ) )              
                ),
        ui.nav_panel("Custom BED",
                    ui.layout_column_wrap(
                        # Input Card
                        ui.card(
                            ui.card_header("Input Panel:"),
                                # Input Files
                                ui.output_ui('update_samples_cb'),
                                ui.output_ui('update_files_cb'),
                                ui.input_file('bed_file_cb', 'Upload BED File'),
                                #Select Method
                                ui.input_select(
                                    'method_cb',
                                    'Choose Depth Extraction Method',
                                    {
                                        'run_pysam_coverage': 'Pysam Coverage',
                                        'run_pysam_depth': 'Pysam depth',
                                        'run_samtools_depth': 'Samtools depth'
                                    }
                                ),
                                ui.input_switch("stats_cb", "Statistics Metrics", False),
                                # Button to run backend scripts to calculate metrics
                                ui.input_action_button('run_cb', label='Run Metrics Analysis'),
                                # Download output
                                ui.input_text('output_cb', 'Output file', placeholder='e.g. metrics.xlsx'),
                                ui.download_button('download_data_cb', 'Download Report'),
                            height='auto'
                        ),
                        # Metrics Display Card
                        ui.card(
                            ui.card_header('Metrics Display:'),
                            popover_style,
                            ui.popover(
                                ui.span(icon, 'Metrics Selection', style='position:absolute; top: 10px; right: 7px; font-size:12px;'),
                                ui.output_ui('metrics_select_cb')
                            ),
                            ui.output_data_frame('display_output_cb'),
                            height='auto'
                        )
                    )
                )
            )

    # About Page
    page_ABOUT = ui.card(
        ui.card_header("About This App"),
        ui.p("This application calculates key quality control (QC) metrics for sequencing alignment files (BAM/CRAM), filtered by user-defined genomic regions provided in a BED file. It helps researchers and bioinformaticians assess data quality within specific regions of interest."),
        ui.p("Version: NGS Metrics Calculator v1.0"),
        ui.p("Developed by: Ana Rita Moreira da Silva, during the MSc internship at Unilabs Genetics."),
        ui.p(
            "Contributors: ",
            ui.tags.a("Pedro Venâncio", href="https://github.com/pfcv99", target="_blank"),
            ", ",
            ui.tags.a("Alberto Pessoa", href="https://github.com/pessoa-am", target="_blank"),
            ", Unilabs Genetics Team"
        ),
        ui.p(
            "Supervised by: ",
            ui.tags.a("Alberto Pessoa", href="https://github.com/pessoa-am", target="_blank"),
            ", Unilabs Genetics Team"
        ),
        ui.p("This work is part of the Master's thesis titled 'Development of quality management software in Next Generation Sequencing (NGS) analysis', completed within the Master in Clinical Bioinformatics (Genome specialization) at the University of Aveiro."),
        ui.p(
            "Contact: ",
            ui.tags.a("anar.m.silva10@gmail.com", href="mailto:anar.m.silva10@gmail.com")
        ),
        ui.h4('Resources'),
        ui.p(
            "If you use this tool in your research or publications, please cite the following resource: ",
            ui.tags.a("GitHub Repository", href="https://github.com/anarmsilva10/metrics_app", target="_blank")
        ),
    )

    return ui.page_fluid(

        ui.panel_title(TITLE) if show_title else None,

        ui.page_navbar(
            ui.nav_panel("Somatics", page_somatics),
            ui.nav_panel("WES", page_WES),
            ui.nav_panel("About", page_ABOUT)
        ),
        theme=shinyswatch.get_theme('united'),
    )

@module.server
def metrics_app_server(input, output, session):

    progression = ui.tags.img(src="https://usagif.com/wp-content/uploads/loading-58.gif",
                              style="width:40px; vertical-align:middle; margin-right:8px;")

    def default_metrics():

        default_selection = ['Average Read Depth', 'Depth of Coverage 1x', 'Depth of Coverage 10x',
                            'Depth of Coverage 15x', 'Depth of Coverage 20x', 'Depth of Coverage 30x', 'Depth of Coverage 50x',
                            'Depth of Coverage 100x', 'Depth of Coverage 500x']
       
        return default_selection

    def get_inputs(input, tab, subtab):

        if tab == "Somatics":

            if subtab == "Somatics":

                req(input.run(), input.cram_file(), input.bed_file(), input.method())

                cram_path = input.cram_file()[0]['datapath']
                bed_path = input.bed_file()[0]['datapath']
                bed_df = filter_bed(bed_path)
                method = input.method()
                stats = input.stats_somatics()

        elif tab == "WES":

            if subtab == "Single Gene":

                # Single Gene
                req(input.run_wes(), input.method_wes())
           
            # why [0], because in s3.py it returns a tuple, where the cram/bam file is the first one
                ui.notification_show(ui.tags.div(progression, "Downloading file..."), id="file", duration=None, type="message")
                cram_path = download_files(input.cram_file_wes())[0]
                bed_df = create_bed(
                    input.genome(),
                    input.gene_select(),
                    input.exon_select_wes(),
                    input.trim_utr_wes(),
                    input.padding_wes()
                )
                method = input.method_wes()
                stats = input.stats_wes()

            # Gene Panel
            elif subtab == "Gene Panel":
                req(input.run_gp(), input.method_gp())
                ui.notification_show(ui.tags.div(progression, "Downloading file..."), id="file", duration=None, type="message")
                cram_path = download_files(input.cram_file_gp())[0]

                if input.gene_select_gp():
                    gene_selection_gp = input.gene_select_gp()
                    if 'exon_select_gp' in input:
                        exon_selection_gp = input.exon_select_gp()
                else:
                    gene_selection_gp = None
                    exon_selection_gp = None

                bed_df = create_panel_bed(
                    input.genome_gp(),
                    gene_selection_gp,
                    exon_selection_gp,
                    input.gene_panels(),
                    input.trim_utr_gp(),
                    input.padding_gp()
                )
                method = input.method_gp()
                stats = input.stats_gp()

            # Custom BED
            elif subtab == "Custom BED":
                req(input.run_cb(), input.bed_file_cb(), input.method_cb())
           
                ui.notification_show(ui.tags.div(progression, "Downloading file..."), id="file", duration=None, type="message")
                cram_path = download_files(input.cram_file_cb())[0]
                bed_path = input.bed_file_cb()[0]['datapath']
                bed_df = filter_bed(bed_path)
                method = input.method_cb()
                stats = input.stats_cb()
       
            ui.notification_remove(id='file')

        else:
            shared_log.logger.warning("Unknown Panel!")
            return None, None, None, None
       
        ref_cache = './data/reference_genome/'
   
        return cram_path, bed_df, method, stats, ref_cache

    #############################
    ### Inputs implementation ###
    #############################

    # Single Gene
    # Sample Dropdown from S3
    @output
    @render.ui
    def update_samples():

        samples = list_samples()

        return ui.input_select('sample_selector', 'Choose Sample:', choices=samples)

    # Files Dropdown from S3
    @output
    @render.ui
    def update_files():

        sample = input.sample_selector()

        if sample:
            files = list_files(sample)

        else:
            files = []

        return ui.input_select('cram_file_wes', 'Choose File:', choices=files)
   
    # Gene Panel
    # Sample Dropdown from S3
    @output
    @render.ui
    def update_samples_gp():

        samples = list_samples()

        return ui.input_select('sample_selector_gp', 'Choose Sample:', choices=samples)

    # Files Dropdown from S3
    @output
    @render.ui
    def update_files_gp():

        sample = input.sample_selector_gp()

        if sample:
            files = list_files(sample)

        else:
            files = []

        return ui.input_select('cram_file_gp', 'Choose File:', choices=files)
   
    # Custom BED
    # Sample Dropdown from S3
    @output
    @render.ui
    def update_samples_cb():

        samples = list_samples()

        return ui.input_select('sample_selector_cb', 'Choose Sample:', choices=samples)

    # Files Dropdown from S3
    @output
    @render.ui
    def update_files_cb():

        sample = input.sample_selector_cb()

        if sample:
            files = list_files(sample)

        else:
            files = []

        return ui.input_select('cram_file_cb', 'Choose File:', choices=files)

    # Single Gene
    # Gene Dropdown
    @reactive.calc
    def gene_choices():
        req(input.genome())

        return get_genes(input.genome())
   
    @reactive.effect
    def gene_selection():
        ui.update_selectize('gene_select', choices=gene_choices())

    # Exon Dropdown for the selected gene
    @output
    @render.ui
    def exon_selector():

        req(input.genome(), input.gene_select())

        exons = get_exons(input.genome(), input.gene_select())

        choices = [exon['number'] for exon in exons]

        return ui.input_selectize('exon_select_wes', 'Select Exons:', choices=choices, multiple=True)
   
    @reactive.effect
    def all_exons():

        req(input.genome(), input.gene_select())

        select_all = input.select_exons()

        selection = get_exons(input.genome(), input.gene_select())

        choices = [exon['number'] for exon in selection]

        ui.update_selectize('exon_select_wes', choices=choices, selected=choices if select_all else [])

    # Gene Panel
    # Gene Panel Files Dropdown
    @output
    @render.ui
    def gene_panel():

        samples = get_gene_panel()

        return ui.input_select('gene_panels', 'Choose Gene Panel:', choices=samples)

    # Gene Dropdown
    @output
    @render.ui
    def select_gene():
        req(input.gene_panels())

        genes = get_panel_genes(input.gene_panels())

        return ui.input_selectize(
            'gene_select_gp',
            'Select Genes:',
            choices=genes,
            multiple=True,
            selected=[]
        )

    # Exon Dropdown for the selected gene
    @output
    @render.ui
    def exon_select():
        req(input.genome_gp(), input.gene_panels(), input.gene_select_gp())

        exons = get_panel_exons(input.genome_gp(), input.gene_select_gp())

        choices = [exon['number'] for exon in exons]

        return ui.input_selectize(
            'exon_select_gp',
            'Select Exons:',
            choices=choices,
            multiple=True
        )

    @reactive.effect
    def all_exons_gp():

        req(input.genome_gp(), input.gene_panels(), input.gene_select_gp())

        select_all = input.select_exons_gp()

        selection = get_panel_exons(input.genome_gp(), input.gene_select_gp())

        choices = [exon['number'] for exon in selection]

        ui.update_selectize('exon_select_gp', choices=choices, selected=choices if select_all else [])

    # Metrics values storage
    metrics_result_somatics = reactive.Value()
    metrics_result_wes = reactive.Value()
    metrics_result_gp = reactive.Value()
    metrics_result_cb = reactive.Value()

    # Parameters storage
    last_params_wes = reactive.Value()
    last_params_gp = reactive.Value()

    ###########################
    ### Metrics Calculation ###
    ###########################

    # For Somatics Page
    @reactive.effect
    @reactive.event(input.run)
    def metrics_calc_somatics():

        cram_path, bed_df, method, stats, ref_cache = get_inputs(input, 'Somatics', 'Somatics')

        if not cram_path or bed_df.empty or not method:
            shared_log.logger.error("Error uploading the required parameters.")
            metrics_result_somatics.set(pd.DataFrame())

        else:
            ui.notification_show(ui.tags.div(progression, "Calculating Metrics..."), id="calc", duration=None, type="message")
            df = calculate_metrics(cram_path, bed_df, None, None, method, stats=stats, ref_cache=ref_cache)
            if stats == True:
                ui.notification_show(ui.tags.div(progression, "Extracting Statistic Metrics..."), id="stats", duration=None, type="message")
            ui.notification_show(ui.tags.div(progression, "Extracting Depth..."), id="depth", duration=None, type="message")
            shared_log.logger.info(f"Metrics Display")
            metrics_result_somatics.set(df)
            ui.notification_remove(id='stats')
            ui.notification_remove(id='depth')
   
        ui.notification_remove(id='calc')

    # For WES - Single Gene Page
    @reactive.effect
    @reactive.event(input.run_wes)
    def metrics_calc_wes():

        cram_path, bed_df, method, stats, ref_cache = get_inputs(input, 'WES', 'Single Gene')

        gene_selection_wes = input.gene_select()
        exon_selection_wes = input.exon_select_wes()

        if not cram_path or bed_df.empty or not method:
            shared_log.logger.error("Error uploading the required parameters.")
            metrics_result_wes.set(pd.DataFrame())
            last_params_wes.set(None)

        else:
            ui.notification_show(ui.tags.div(progression, "Calculating Metrics..."), id="calc", duration=None, type="message")
            df = calculate_metrics(cram_path, bed_df, gene_selection=gene_selection_wes, exon_selection=exon_selection_wes, method=method, stats=stats, ref_cache=ref_cache)
            if stats == True:
                ui.notification_show(ui.tags.div(progression, "Extracting Statistic Metrics..."), id="stats", duration=None, type="message")
            ui.notification_show(ui.tags.div(progression, "Extracting Depth..."), id="depth", duration=None, type="message")
            shared_log.logger.info(f"Metrics Display")
            metrics_result_wes.set(df)
            ui.notification_remove(id='stats')
            ui.notification_remove(id='depth')

            last_params_wes.set({
                'cram_path': cram_path,
                'bed_df': bed_df,
                'method': method
            })
       
        ui.notification_remove(id='calc')

    # For WES - Gene Panel Page
    @reactive.effect
    @reactive.event(input.run_gp)
    def metrics_calc_gp():

        cram_path, bed_df, method, stats, ref_cache = get_inputs(input, 'WES', 'Gene Panel')
       
        if input.gene_select_gp():
            gene_selection_gp = input.gene_select_gp()
            if 'exon_select_gp' in input:
                exon_selection_gp = input.exon_select_gp()
        else:
            gene_selection_gp = None
            exon_selection_gp = None

        if not cram_path or bed_df.empty or not method:
            shared_log.logger.error("Error uploading the required parameters.")
            metrics_result_gp.set(pd.DataFrame())
            last_params_gp.set(None)

        else:
            ui.notification_show(ui.tags.div(progression, "Calculating Metrics..."), id="calc", duration=None, type="message")
            df = calculate_metrics(cram_path, bed_df, gene_selection_gp, exon_selection_gp, method, stats=stats, ref_cache=ref_cache)
            if stats == True:
                ui.notification_show(ui.tags.div(progression, "Extracting Statistic Metrics..."), id="stats", duration=None, type="message")
            ui.notification_show(ui.tags.div(progression, "Extracting Depth..."), id="depth", duration=None, type="message")
            shared_log.logger.info(f"Metrics Display")
            metrics_result_gp.set(df)
            ui.notification_remove(id='stats')
            ui.notification_remove(id='depth')
            last_params_gp.set({
                'cram_path': cram_path,
                'bed_df': bed_df,
                'method': method
            })
       
        ui.notification_remove(id='calc')

    # For WES - Custom BED Page
    @reactive.effect
    @reactive.event(input.run_cb)
    def metrics_calc_cb():

        cram_path, bed_df, method, stats, ref_cache = get_inputs(input, 'WES', 'Custom BED')

        if not cram_path or bed_df.empty or not method:
            shared_log.logger.error("Error uploading the required parameters.")
            metrics_result_cb.set(pd.DataFrame())

        else:
            ui.notification_show(ui.tags.div(progression, "Calculating Metrics..."), id="calc", duration=None, type="message")
            df = calculate_metrics(cram_path, bed_df, None, None, method, stats=stats, ref_cache=ref_cache)
            if stats == True:
                ui.notification_show(ui.tags.div(progression, "Extracting Statistic Metrics..."), id="stats", duration=None, type="message")
            ui.notification_show(ui.tags.div(progression, "Extracting Depth..."), id="depth", duration=None, type="message")
            shared_log.logger.info(f"Metrics Display")
            metrics_result_cb.set(df)
            ui.notification_remove(id='stats')
            ui.notification_remove(id='depth')

        ui.notification_remove(id='calc')

    #########################
    ### Metrics Selection ###
    #########################

    def get_metrics(tab, subtab):

        if tab == "Somatics":

            if subtab == "Somatics":
                metrics_df = metrics_result_somatics.get()

        elif tab == "WES":

            if subtab == "Single Gene":
                metrics_df = metrics_result_wes.get()

            elif subtab == "Gene Panel":
                metrics_df = metrics_result_gp.get()

            elif subtab == "Custom BED":  
                metrics_df = metrics_result_cb.get()

        else:
            shared_log.logger.warning("Unknown Panel!")
            return None, None
       
        if metrics_df is None or metrics_df.empty:
            choices = []
        else:
            choices = list(metrics_df.columns)

        default_selection = [m for m in default_metrics() if m in choices]

        return choices, default_selection
   
    # Somatics Page
    @output
    @render.ui
    def metrics_select_somatics():

        choices, default_selection = get_metrics('Somatics', 'Somatics')

        return ui.div(
            ui.input_checkbox_group(
            'metrics_selection_somatics',
            'Select metrics to display:',
            choices = choices,
            selected = default_selection
        ), class_='metrics')
   
    # WES - Single Gene Page
    @output
    @render.ui
    def metrics_select_wes():

        choices, default_selection = get_metrics('WES', 'Single Gene')

        return ui.div(
            ui.input_checkbox_group(
            'metrics_selection_wes',
            'Select metrics to display:',
            choices = choices,
            selected = default_selection
        ), class_='metrics')
   
    # WES - Gene Panel Page
    @output
    @render.ui
    def metrics_select_gp():

        choices, default_selection = get_metrics('WES', 'Gene Panel')

        return ui.div(
            ui.input_checkbox_group(
            'metrics_selection_gp',
            'Select metrics to display:',
            choices = choices,
            selected = default_selection
        ), class_='metrics')
   
    # WES - Custom BED Page
    @output
    @render.ui
    def metrics_select_cb():

        choices, default_selection = get_metrics('WES', 'Custom BED')

        return ui.div(
            ui.input_checkbox_group(
            'metrics_selection_cb',
            'Select metrics to display:',
            choices = choices,
            selected = default_selection
        ), class_='metrics')

    #######################
    ### Metrics Display ###
    #######################

    def get_dataframe(tab, subtab):
       
        if tab == "Somatics":

            if subtab == "Somatics":
                metrics_df = metrics_result_somatics.get()

        elif tab == "WES":

            if subtab == "Single Gene":
                metrics_df = metrics_result_wes.get()

            elif subtab == "Gene Panel":
                metrics_df = metrics_result_gp.get()

            elif subtab == "Custom BED":
                metrics_df = metrics_result_cb.get()
        else:
            shared_log.logger.warning("Unknown Panel!")
            return pd.DataFrame({"Message": ["Unknown panel selected."]})

        if metrics_df is None or metrics_df.empty:
            return pd.DataFrame({"Message": ["No Data. Select BAM and BED Files."]})

        try:
            selected_metrics = select_metrics(tab, subtab)
        except Exception:
            selected_metrics = None

        if not selected_metrics:
            selected_metrics = [m for m in default_metrics() if m in metrics_df.columns]

        selected_metrics = list(selected_metrics)

        try:
            chose_df = metrics_df[selected_metrics]
        except KeyError:
            return pd.DataFrame({"Message": ["Selected metrics not found in dataframe."]})
       
        transposed_df = chose_df.T

        if 'Gene' in metrics_df.columns or 'Exon' in metrics_df.columns:
            # Assign column names from the first column of metrics_df
            new_cols = metrics_df.iloc[:, 0].astype(str).tolist()
            if len(new_cols) == len(transposed_df.columns):
                transposed_df.columns = new_cols
        else:
            transposed_df.columns = ['Results']

        transposed_df = transposed_df.reset_index().rename(columns={'index': 'Metric'})
        transposed_df.columns = transposed_df.columns.map(str)

        cols_to_convert = transposed_df.columns.difference(['Metric'])
        transposed_df[cols_to_convert] = transposed_df[cols_to_convert].round(0)

        return transposed_df
   
    def select_metrics(tab, subtab):

        if tab == "Somatics":

            if subtab == "Somatics":        
                selected_metrics = input.metrics_selection_somatics()

        elif tab == "WES":

            if subtab == "Single Gene":      
                selected_metrics = input.metrics_selection_wes()

            elif subtab == "Gene Panel":        
                selected_metrics = input.metrics_selection_gp()

            elif subtab == "Custom BED":        
                selected_metrics = input.metrics_selection_cb()

        else:
            shared_log.logger.warning("Unknown Panel!")
            return None

        if not selected_metrics:
            return None

        return list(selected_metrics)

    # Somatics Page Output
    @output
    @render.data_frame
    def display_output_somatics():

        transposed_df = get_dataframe('Somatics', 'Somatics')

        return render.DataGrid(transposed_df, width='1000px')
   
    # WES - Single Gene Page Output
    @output
    @render.data_frame
    def display_output_wes():

        transposed_df = get_dataframe('WES', 'Single Gene')

        return render.DataGrid(transposed_df, width='1000px')
   
    # WES - Gene Panel Page Output
    @output
    @render.data_frame
    def display_output_gp():

        transposed_df = get_dataframe('WES', 'Gene Panel')

        return render.DataGrid(transposed_df, width='1000px')
   
    # WES - Custom BED Page Output
    @output
    @render.data_frame
    def display_output_cb():

        transposed_df = get_dataframe('WES', 'Custom BED')

        return render.DataGrid(transposed_df, width='1000px')
   
    ########################
    ### Metrics Download ###
    ########################

    def writer(tab, subtab):

        shared_log.logger.info("Downloading Metrics Report")

        if tab == "Somatics":

            if subtab == "Somatics":
                metrics_df = metrics_result_somatics.get()

        elif tab == "WES":

            if subtab == "Single Gene":
                metrics_df = metrics_result_wes.get()

            # Gene Panel
            elif subtab == "Gene Panel":
                metrics_df = metrics_result_gp.get()

            # Custom BED
            elif subtab == "Custom BED":
                metrics_df = metrics_result_cb.get()

        else:
            shared_log.logger.warning("Unknown Panel!")
            return None, None
       
        if metrics_df.empty:
            return
           
        output = io.BytesIO()

        metrics_df.to_excel(output, index=False)
        output.seek(0)

        yield output.read()

    # Somatics Page - download
    @output
    @render.download(filename=lambda: input.output_somatics())
    def download_data_somatics():

        return writer('Somatics', 'Somatics')

    # WES - Single Page download
    @output
    @render.download(filename=lambda: input.output_wes())
    def download_data_wes():

        return writer('WES', 'Single Gene')
   
    # WES - Gene Panel Page download
    @output
    @render.download(filename=lambda: input.output_gp())
    def download_data_gp():

        return writer('WES', 'Gene Panel')
   
    # WES - Custom BED Page download
    @output
    @render.download(filename=lambda: input.output_cb())
    def download_data_cb():

        return writer('WES', 'Custom BED')
   
    #####################
    ### Coverage Plot ###
    #####################
   
    # WES - Single Gene
    @output
    @render_widget
    def plot_coverage_wes():
        if not input.plot_wes():
            return go.Figure()
       
        params = last_params_wes.get()
        if not params:
            return go.Figure()
       
        bed_df = params['bed_df'].copy()
       
        selected_exons = input.plot_exons_wes()
        if selected_exons:
            bed_df = bed_df[bed_df['EXON'].isin(selected_exons)]

        fig = coverage_plot(
            cram_path=params['cram_path'],
            bed_df=bed_df,
            method=params['method'],
            threshold=input.threshold_wes(),
            padding=input.padding_wes()
        )
       
        return fig

    @reactive.effect
    def update_plot_exons_choices():

        req(input.gene_select())

        exons = get_exons(input.genome(), input.gene_select())
        choices = [exon['number'] for exon in exons]

        ui.update_selectize('plot_exons_wes', choices=choices, selected=[])

    # WES - Gene Panel
    @output
    @render_widget
    def plot_coverage_gp():
        if not input.plot_gp():
            return go.Figure()

        params = last_params_gp.get()
        if not params:
            return go.Figure()

        bed_df = params['bed_df'].copy()
       
        gene = input.gene_select_plot_gp()
        bed_df = bed_df[bed_df['GENE'] == gene]
       
        selected_exons = None
        if hasattr(input, 'plot_exons_gp'):
            try:
                selected_exons = input.plot_exons_gp()
            except Exception:
                selected_exons = None

        if selected_exons:
            bed_df = bed_df[bed_df['EXON'].isin(selected_exons)]

        fig = coverage_plot(
            cram_path=params['cram_path'],
            bed_df=bed_df,
            method=params['method'],
            threshold=input.threshold_gp(),
            padding=input.padding_gp()
        )
        return fig

    @output
    @render.ui
    def gene_select_plotly_gp():
        req(input.gene_panels())

        selected_genes_metrics = input.gene_select_gp()
        if not selected_genes_metrics:
            genes = get_panel_genes(input.gene_panels())

            return ui.input_selectize(
                'gene_select_plot_gp',
                'Select Gene to Plot:',
                choices=genes,
                multiple=False,
                selected=None
            )
        else:
           
            return ui.input_selectize(
                'gene_select_plot_gp',
                'Select Gene to Plot:',
                choices=selected_genes_metrics,
                multiple=False,
                selected=selected_genes_metrics[0]
            )

    @output
    @render.ui
    def plot_exons_selector_gp():
        req(input.gene_select_plot_gp())

        selected_exon_metrics = None
        if hasattr(input, 'exon_select_gp'):
            try:
                selected_exon_metrics = input.exon_select_gp()
            except Exception:
                selected_exon_metrics = None

        if not selected_exon_metrics:
            exons = get_panel_exons(input.genome_gp(), [input.gene_select_plot_gp()])
            choices = [exon['number'] for exon in exons]
            return ui.input_selectize(
            'plot_exons_gp',
            'Select Exons to Plot:',
            choices=choices,
            selected=[],
            multiple=True
            )
        else:
            choices = selected_exon_metrics
            return ui.input_selectize(
            'plot_exons_gp',
            'Select Exons to Plot:',
            choices=choices,
            selected=choices[0],
            multiple=True
            )

atexit.register(clean_temp_dirs)

app = App(
    ui=metrics_app_ui(NS, show_title=True),
    server=lambda input, output, session: metrics_app_server(NS)
)