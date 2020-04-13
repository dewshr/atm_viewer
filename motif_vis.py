#######################################################################
#Author: Dewan Shrestha                                               #
#PI    : Yong Cheng (yong.cheng@stjude.org)                           #              
#Email : dewshrs@gmail.com / dshrest2@uthsc.edu / dshresth@stjude.org #
#######################################################################

import sys
import os
script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_path)
sys.path.append(script_path+'/scripts')
from utils import *
import argparse
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import webbrowser
from dash.exceptions import PreventUpdate


############################################### initial arguments ##############################
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default=None, help='multiple alignment format file for visualization')
parser.add_argument('-d','--dir', default='motif_vis_temp_files', help='folder name for the temporary files generated during the process')
parser.add_argument('-l','--load', default = True, help='automatically opens browser if value is true')
parser.add_argument('-f','--format', default='maf', choices=['maf','clustal', 'fasta'],help='format of the alignment file')

args = parser.parse_args()

#creating log file
logger_path = os.path.abspath(args.dir)
logger.add(logger_path+'/visualization_{time}.log', rotation='10 MB')

#if output directory is not present, creates a new directory
if not os.path.exists(args.dir):
	os.makedirs(args.dir)


if args.input == None:
	logger.info('using default example')
	extract_seq_from_maf(script_path+'/example/example.maf', args.dir)
	aln = AlignIO.read(os.path.join(args.dir,'temp_aln.fa'),'fasta')
	#parser.print_help()
	#sys.exit()
else:
	if args.format=='maf':
		extract_seq_from_maf(args.input, args.dir) # processing the maf file to combine alignment blocks
		aln = AlignIO.read(os.path.join(args.dir,'temp_aln.fa'),'fasta')
	elif args.format == 'clustal':
		extract_seq_from_clustal(args.input, args.dir)
		aln = AlignIO.read(os.path.join(args.dir,'temp_aln.fa'),'fasta')
	elif args.format == 'fasta':
		aln = AlignIO.read(args.input,'fasta')
		with open(os.path.join(args.dir ,'fimo_input.fa'),'w') as f:
			file_ = open(args.input).read().replace('-','')
			f.write(file_)
	else:
		logger.info('the format is not supported, please use the supported format file')
		sys.exit()


			

#creating lists for tf and species for dash dropdown options
motif_options = [{'label':'{} ({})'.format(jasper_dict[key]['motif_name'],jasper_dict[key]['species'] ), 'value':key} for key,value in jasper_dict.items()]
species_options= [{'label':value, 'value':key} for key, value in species_dict.items()]
logger.info('options created for the buttons')



#extracting sequences and ids from the alignment file
#aln = AlignIO.read(os.path.join(args.dir,'temp_aln.fa'),'fasta') 
seqs, seq_ids, nucleotide_bases, seq_names, sp_names  = [], [], [], [], []

for rec in aln:
	seqs.append(rec.seq)
	seq_ids.append(rec.id.lower())


for s in seq_ids:
	for sp in species_options:
		if sp['value'].lower() in s:
			seq_names.append([sp['label'], s])
			sp_names.append(sp['label'])

logger.info('all the species found')
logger.info(sp_names)
for s in seqs:
	nucleotide_bases.append([i for i in s])
#nucleotide_bases2 = [item for sublist in nucleotide_bases for item in sublist]
seq_len = len(nucleotide_bases[0])
default_base_values = [[0]*seq_len for x in range(len(nucleotide_bases))]
hover_data= [['']*seq_len for x in range(len(nucleotide_bases))]

#y = [[i[0]]*seq_len for i in seq_names]
#y = [item for sublist in y for item in sublist]






app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])


################################## top header ########################################
top_header= dbc.Row(dbc.Col(html.H1('Motif Turnover Visualization',
							style={'textAlign':'center',
								   'marginTop':5,
								   'marginBottom':10,
								   'paddingBottom':5,
								   'fontFamily':'Courier New',
								   'fontWeight':'bold',
								   'fontColor':'#FF0000'
								  })))




################################### empty row at the end #############################
empty_row = dbc.Row(dbc.Col(html.P('.',style={'marginTop':10,
											 'fontSize':5})))



################################### dropdown options #############################################
dropdown = html.Div([
					html.H3('Select the TF you want to search for:', 
							style={'marginTop':15,
									'paddingLeft':10,
								   'fontSize':20,
								   'fontFamily':'Courier New',
								   'fontWeight':'bold'
								  }),
	
					 dcc.Dropdown(id='motif_id',
								options=motif_options,
								multi=True,
								value=['MA0139.1'],
								 style={'marginTop':5,
									   'paddingTop':5,
									   'paddingBottom':5,
									   'paddingLeft':10}
								 
								),
	
					html.Button(id='submit', n_clicks=0, children='Submit', 
								style={'fontSize':20,
									  'paddingLeft':10,
									   'marginLeft':10,
									  'fontFamily':'Courier New'}),
	
					dcc.Loading(id='fimo_running',
							   children=[html.Div(id='fimo_status',
												  style={'paddingLeft':10,
														 'paddingTop':10
														 
															   
							   })]),
	
					html.H3('Select the species you want to view:',
							style={'marginTop':20,
								   'marginBottom':5,
								   'fontSize':20,
								   'fontFamily':'Courier New',
								   'fontWeight':'bold',
								   'paddingLeft':10
								  }),
	
					dcc.Dropdown(id='species_id',
								options=species_options,
								multi=True,
								value=['hg','mm', 'pantro', 'gorgor', 'rhemac'],
								 style={'marginTop':10,
									   'paddingTop':5,
									   'paddingBottom':5,
									   'paddingLeft':10}
								)
							
				   
					
					])




######################################### figure bar legend #######################################
heatmap_legend = html.Div([
							dcc.Graph(id='figure_bar',
										config={'toImageButtonOptions':{'format':'svg'},
											   'displayModeBar':True}
										 )
							], style={'paddingRight':10,
										'paddingLeft':0,
										'marginLeft':0
									})
	

tree = html.Div([
			dcc.Graph(id='tree',
						config={'displayModeBar':False}
								)
	],style={
			'paddingRight':0,
			'marginRight':0
	})

###################################### [motif informations: seleected vs found] and [species] ############
tf_selected =  [
	dbc.CardHeader("TF Selected"),
	dbc.CardBody(
		[
			#html.H5("Card title", className="card-title"),
			dcc.Markdown(id= 'tf_selected',
				className="card-text",
			)
		]
	)
]
tf_found = [
	dbc.CardHeader("TF Found"),
	dbc.CardBody(
		[
			#html.H5("Card title 2", className="card-title"),
			dcc.Markdown(id='tf_found_ids',
				className="card-text",
			),
		]
	),
]

species_details = [
	dbc.CardHeader("Details about species aligned"),
	dbc.CardBody(
		[
			#html.H5("Card title 3", className="card-title"),
			dcc.Markdown(id='species',children = get_markdown(seq_names, fimo=False),
				className="card-text",
			),
		]
	),
]



################################## cards with information about tfs found and species #######################
cards = html.Div(
	[
		dbc.Row(
			[
				dbc.Col(dbc.Card(tf_selected, color="info", inverse=True)),
				dbc.Col(dbc.Card(tf_found, color="success", inverse=True))
			],
			className="mb-4"
		),
		dbc.Row(
			[
				dbc.Col(dbc.Card(species_details, color="secondary", inverse=True)),
			],
			className="mb-4")
	], style={
		'paddingLeft':10
	})

description = html.Div(cards)






################################################# main alignment file ###############################
main_figure = html.Div([
						dcc.Graph(id='motif_plot',
									config={'toImageButtonOptions':{'format':'svg'},
											'displayModeBar':True
											}
										 )
										],style={
										'paddingRight':10,
										'marginBottom':30
											})

end =dbc.Row(dbc.Col(html.P('.',style={'marginTop':10,
											 'fontSize':5})))





################################################### layout ########################################
body = html.Div(dbc.Row([
				dbc.Col(dropdown, width=4),
				dbc.Col(tree, width=3),
				dbc.Col(heatmap_legend, width=5)
]))

graph = html.Div(dbc.Row([
				dbc.Col(description, width=4),
				dbc.Col(main_figure, width=8)
]))

#dcc.Store(id='store_data')
app.layout=html.Div([dcc.Store(id='memory-output'),
			top_header, empty_row,body, graph, end
]
)


######################################################### callbacks #################################
@app.callback([Output('tf_selected','children'),
			   Output('tf_found_ids','children'),
			   Output('memory-output', 'data'),
			   Output('fimo_status','children')
			  ],
			 [Input('submit','n_clicks')],
			 [State('motif_id','value')])
def tf_selected(n_clicks, value):
	tf_selected = get_markdown(value)
	if n_clicks > 0:
		logger.info('\n\n------------number of clicks{}------------\n\n', n_clicks)
		alert =  dbc.Alert("Fimo run completed", color="success")
		logger.info('running fimo')
		pwm_file = extract_pwm(value, args.dir)
		os_code, fimo_file = run_fimo(os.path.join(args.dir,'fimo_input.fa'), pwm_file, args.dir)

		if os_code == 1:
			logger.info('Something went wrong while running fimo')
			sys.exit(1)
		else:
			logger.info('processing fimo')
			fimo_ids, fimo_list, motif_species = process_fimo(fimo_file)
			
			logger.info('fimo_ids\n{}',fimo_ids)
			logger.info('fimolist\n{}', fimo_list)
			
			tf_found_ = get_markdown(fimo_ids)
			logger.info('\n\n------------ NUMBER OF TF FOUND : {}------------\n\n', len(fimo_ids))
			logger.info('tf found\n{}', tf_found_)
			#print(hover_data[0])
			base_color_values, hover_values, colorscale, legend_col = get_intermediate_data(fimo_ids, fimo_list, default_base_values, hover_data, nucleotide_bases)
			#logger.info('hover values \n{}',hover_values)
			#logger.info('\nhover_data \n {}', hover_data)
			#logger.info('colorscale\n{}',colorscale)
			#logger.info('legend_col\n{}', legend_col)
			data = [{'base_values':base_color_values}, {'hover_values':hover_values},{'nucleotide_bases':nucleotide_bases},{'sp_names':sp_names}, {'colorscale':colorscale},{'legend_col':legend_col}, {'motif_species':motif_species}]
	else:
		tf_found_ =[None]
		alert = dbc.Alert("Press submit to run fimo", color="info")
		data=None
	#'''.format(motif_names[value[0]], motif_names[value[1]])
	#return 'The motifs selected are: {}'.format(', '.join([jasper_dict[i]['motif_name'] for i in value]))
	return tf_selected, tf_found_ , data, alert

@app.callback([Output('motif_plot','figure'),
			   Output('figure_bar','figure'),
			   Output('tree', 'figure')
			  ],
			 [Input('species_id','value'),
			 Input('memory-output','data')])
def create_motif_plot(species_selected, data):
	if data is None:
		raise PreventUpdate

	indices=[]
	#print('species_selected',species_selected)
	for n in species_selected:
		try:
			indices.append(sp_names.index(species_dict[n]))
		except:
			continue
	#print(indices)
	return get_figure(data, indices)

if __name__ == '__main__':
	if args.load == True:
		webbrowser.open_new('http://127.0.0.1:8051')
		app.run_server(debug=True, port=8051)
	else:
		app.run_server(debug=True, port=8051)


















