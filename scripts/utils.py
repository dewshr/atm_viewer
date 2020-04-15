######################################################################
#Author: Dewan Shrestha                                              #
#PI    : Yong Cheng (yong.cheng@stjude.org)                          #
#Email: dewshrs@gmail.com / dshrest2@uthsc.edu / dshresth@stjude.org #
######################################################################

import numpy as np
import os
import sys
import plotly
import plotly.graph_objs as go
import json
import re
from Bio import AlignIO
from Bio.Seq import Seq
from loguru import logger



data_path = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1])
#print(script_path)
# list of color assuming the user will select maximum of 14 tfs at time
color_list=['whitesmoke','limegreen','tomato','teal','lightpink','steelblue','darkmagenta','orange','mediumslateblue','darkred','forestgreen','goldenrod','mediumblue','dimgray','darkturquoise']


#header of pwm matrix for meme format
header ='MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\n\nA 0.25 C 0.25 G 0.25 T 0.25\n\n'



#reading in dictionary of tf_id, tf_name and pwm
logger.info('loading jasper_dict from jasper_2020_core_vertebrates_insects_non-redundant_pfms_meme.json')
with open(data_path+'/data/jasper_2020_core_vertebrates_insects_non-redundant_pfms_meme.json') as f:
	jasper_dict = json.load(f)



#reading in json file which contains species annotation
logger.info('loading dictionary for species annotation')
with open(data_path+'/data/species_annotation.json') as f:
	species_dict = json.load(f)


''' 
this function takes the multiple alignment format file which contains the alignments in blocks for different organisms 
and generates files with combined sequences from different blocks. It generates two files: one alignment (temp_aln.fa) file 
with dashes and another without dashes (fimo_input.fa) that is used to run fimo to find the motifs.
'''
@logger.catch
def extract_seq_from_maf(maf_file, dir):
	maf_file_data = open(maf_file)
	maf_aln={}
	for line in maf_file_data:
		if line.startswith('s') and 'scaffold' not in line:
			data = line.split()
			species = data[1]
			seq= data[-1].replace('\n','')
			length= len(seq)
			start = int(data[2])
			stop = start + int(data[3])
			strand = data[4]
			#print(start)
			if species not in maf_aln.keys():
				maf_aln[species] ={'start':[start], 'stop':[stop],'strand':[strand],'seq':seq}
				#print(maf_aln)
				
			maf_aln[species]['start'] = maf_aln[species]['start'] + [start]
			maf_aln[species]['stop'] = maf_aln[species]['stop'] + [stop]
			maf_aln[species]['strand'] = maf_aln[species]['strand'] + [strand]
			maf_aln[species]['seq'] = maf_aln[species]['seq']+ seq


	logger.info('creating input file for running fimo: fimo_input.fa and alignment file: temp_aln.fa for visualization')
	fimo_input = open(os.path.join(dir, 'fimo_input.fa'),'w')
	with open(os.path.join(dir ,'temp_aln.fa'),'w') as f:
		for key, value in maf_aln.items():
			f.write('>{}.{}:{}.{}\n{}\n'.format(key, min(value['start']), max(value['stop']), set(value['strand']), value['seq']))
			fimo_input.write('>{}.{}:{}.{}\n{}\n'.format(key, min(value['start']), max(value['stop']), set(value['strand']), value['seq'].replace('-','')))
	
	fimo_input.close()

################### this function will be called if the file is in clustal alignment format ##########
@logger.catch
def extract_seq_from_clustal(clustal_file, dir):
	align = AlignIO.read(clustal_file, 'clustal')
	fimo_input = open(os.path.join(dir, 'fimo_input.fa'),'w')
	with open(os.path.join(dir ,'temp_aln.fa'),'w') as f:
		for rec in align:
			fimo_input.write('>{}\n{}\n'.format(str(rec.id), str(rec.seq).replace('-','')))
			f.write('>{}\n{}\n'.format(str(rec.id),str(rec.seq)))
	fimo_input.close()



############################## extract pwm based on given list of tfs ############################
@logger.catch
def extract_pwm(tf_list, dir):
	with open(os.path.join(dir,'tf_selected_pwm.txt'), 'w') as f:
		f.write(header)
		for tf in tf_list:
			#print(jasper_dict[tf])
			f.write(jasper_dict[tf]['pwm'])

	return os.path.join(dir,'tf_selected_pwm.txt')





############################### run fimo for the selected TFs ##########################
@logger.catch
def run_fimo(file, pwm, dir):

	if os.system('fimo -oc {} {} {}'.format(os.path.join(dir, 'fimo_output'), pwm, file)) ==0:
		logger.info('running fimo completed')
		return 0, os.path.join(dir, 'fimo_output/fimo.tsv')
	else:
		logger.info('Problem running fimo. Please check meme is installed properly')
		sys.exit()
		return 1, os.path.join(dir, 'fimo_output/fimo.tsv')




############################## Process the fimo output to generate the list of found TFs: fimo_ids and motif sequence associated with it (fimo_list)###########		

@logger.catch
def process_fimo(fimo_file):
	logger.info('generating fimoids and fimolist')
	fimo_ids = []
	fimo_list = []
	motif_species= {} #stores TF and species which has that TF
	file_ = open(fimo_file)
	
	for line in file_:
		if line.startswith('M'):
			data = line.split()
			if data[1] not in motif_species.keys():
				motif_species[data[1]] = data[2]
			motif_species[data[1]] = motif_species[data[1]] + data[2]
			fimo_ids.append(data[0])
			if data[5] == '-':
				motif_seq = str(Seq(data[-1]).reverse_complement())
			else:
				motif_seq = data[-1]
			fimo_list.append([motif_seq, data[0]])
			
	fimo_ids = list(set(fimo_ids))
	fimo_ids.sort()
	fimo_list = [list(item) for item in set(tuple(row) for row in fimo_list)]
	fimo_list.sort( key=lambda l: (len(l[0]),l[1]), reverse=True)

	return fimo_ids, fimo_list, motif_species





######################## assigning the color values to the positions where the motif sequence is found
@logger.catch
def motif_color(nucleotide_bases,motif_details, base_values, hover_data):
	hover_values = [i.copy() for i in hover_data.copy()]
	base_color_values = [i.copy() for i in base_values.copy()]
	for n in range(len(nucleotide_bases)):
		print(n)
		color_tracker = {}
		for motif,val in motif_details.items():
			iter = re.finditer(r"{}".format('-*'.join(list(motif))), ''.join(nucleotide_bases[n])) # searching for the motif
			indices = [[m.start(0),m.end(0)] for m in iter] # stores indices for the searched motif if found
			#logger.info(val['mname'])
			#logger.info(indices)
			for index in indices:
				for i in range(index[0], index[1]):
					#logger.info(i)
					#logger.info('changing values now')
					#base_color_values[n][i]= val['mval'] #assigning the value associated with respective TF
					#hover_values[n][i] = val['mname']
					if hover_values[n][i] == '':
						#logger.info(hover_values[n][i])
						base_color_values[n][i]= val['mval']
						hover_values[n][i] = val['mname'] # assigning the motif name for annotation in plot
						color_tracker[i] =1
						#logger.info(hover_values[n][i])
					else:
						r = color_tracker[i] * 0.04
						if color_tracker[i] ==3:
							r = (2 * 0.04) + 0.02
						color_tracker[i] = color_tracker[i] + 1
						if val['mval'] > base_color_values[n][i]:
							base_color_values[n][i]= val['mval'] - r
						else:
							base_color_values[n][i]= base_color_values[n][i] - r
							
						hover_values[n][i] = hover_values[n][i] + ' : '+val['mname'] # combining the name if two motifs are found on same position
	#logger.info(base_color_values)
	#logger.info(hover_values[0])
	return base_color_values, hover_values





####################################### motif dictionary with details regarding color, id, seq and name #############################
@logger.catch
def get_motif_details_dict(fimo_list, legend_col_val, id_color):
	motif_details ={}
	for x in fimo_list:
		if x[0] not in motif_details.keys():
			motif_details[x[0]]= {'mid':x[1],'mval':legend_col_val[x[1]], 'mname':jasper_dict[x[1]]['motif_name'], 'mcolor':id_color[x[1]]}
		else:
			motif_details[x[0]]['mname'] = motif_details[x[0]]['mname'] +' : '+ jasper_dict[x[1]]['motif_name']
	return motif_details





############################# generating data for seuquence color(base_values), sequence annotation (hover_values) ###########################
@logger.catch
def get_intermediate_data(fimo_ids, fimo_list, base_values, hover_data, nucleotide_bases):
	colors_ = color_list[1:len(fimo_ids)+1]
	id_color = {fimo_ids[i]:colors_[i] for i in range(len(fimo_ids))}
	
	if len(fimo_ids) ==0:
		colorscale =[[0,'whitesmoke']]
		legend_col_val={}
		legend_col=[]

		return base_values, hover_data, colorscale, legend_col
	else:
		val = 1.0/(len(fimo_ids))
		if len(fimo_ids) ==1:
			val_list=[1.0]
		elif len(fimo_ids) == 2:
			val_list= [0.45,1]
		else:
			val_list = np.arange(val,1.01,val).round(decimals=2)

		colorscale = [[0,'whitesmoke']]+[list(x) for x in zip(val_list, colors_)]
		logger.info('colorscale \n{}',colorscale)
		legend_col_val = {fimo_ids[i]:val_list[i] for i in range(len(fimo_ids))}
		logger.info('\nlegend_col_val\n {}', legend_col_val)
		
		#motif_details={x[0]:{'mid':x[1],'mval':legend_col_val[x[1]], 'mname':jasper_dict[x[1]]['motif_name'], 'mcolor':id_color[x[1]]} for x in fimo_list}
		motif_details = get_motif_details_dict(fimo_list, legend_col_val, id_color)
		logger.info('motif_details\n{}', motif_details)

		legend_col= [[val['mname'],val['mcolor']]for key,val in motif_details.items()]
		legend_col = list(set(tuple(x) for x in legend_col))

		logger.info('legend_col\n{}',legend_col)
		#print('\n--------------------\n')
		#print(colorscale)
		
		base_color_values, hover_values = motif_color(nucleotide_bases,motif_details, base_values, hover_data)
		

		return base_color_values, hover_values, colorscale, legend_col

	



###################### creates markdown format of the given list ##############################
@logger.catch
def get_markdown(list_, fimo=True):
	text=''
	if fimo==True:
		for item in list_:
			text = text+ '\n- '+ jasper_dict[item]['motif_name']
			#print(text)
		return text
	else:
		for names in list_: # fimo = False, means its species listt
			text = text + '\n- {} :\t{}\n'.format(names[0], ' '.join(names[1].split('.'))) 
		return text



####################### convert scientific name to common name for species ##################
@logger.catch
def get_key(val): 
    for key, value in species_dict.items(): 
         if val == value: 
             return key 	
   
	

############################# getting overall view for the found TF and selected species #####################
#sp_names seq_names_ get_key(value)
@logger.catch
def tree_diagram(seq_names_, motif_species, sp__):
	t=0
	#sp__ = [get_key(x) for x in seq_names_]
	lp = len(sp__)
	tree = go.Figure()
	annotations=[]
	for key, value in motif_species.items():
		x=[t]*lp
		y_val= np.arange(0, lp, 1)
		color=['tomato' if x in value.lower() else 'grey' for x in sp__]
		a1 = dict(x=t,y=y_val[-1]+1,xref='x',yref='y',
            text= key,#'<font size="10"><b>'+key+'</b></font>',
            showarrow=False,
            align='center',
            textangle=-90,
            font=dict(family='Courier New',
            size=15)
            )
		annotations.append(a1)
		tree.add_trace(go.Scatter(x=x, y=y_val, 
                        mode='markers',
                        hoverinfo='none',
                        marker=dict(size=15,
                                color=color   
                                ))
					)

		t=t+1
	#return tree

	tree.update_layout(
            autosize=True,
            showlegend=False,
            annotations=annotations,
            #width=len(motif_species)*100+100,
            #height= (lp__*60) +50,
            plot_bgcolor='white',
            yaxis={'ticktext':seq_names_,'tickvals':y_val},
            xaxis ={'showticklabels':False}
			)

	return tree


######################################### generating alignment plot and figure bar with legend ########################
@logger.catch
def get_figure(data, indices):
	#[{'base_values':base_values}, {'hover_values':hover_values},{'nucleotide_bases':nucleotide_bases},{'sp_names':sp_names}, {'colorscale':colorscale},{'legend_col':legend_col}]
	indices= indices[::-1]
	fig_base_values = data[0]['base_values']
	fig_hover_values = data[1]['hover_values']
	fig_nucleotide_bases = data[2]['nucleotide_bases']
	sp_names = data[3]['sp_names']
	colorscale = data[4]['colorscale']
	legend_col = data[5]['legend_col']
	motif_species = data[6]['motif_species']

	seq_len = len(fig_nucleotide_bases[0])


	logger.info(motif_species)
	logger.info(legend_col)
	logger.info('colorscale')
	logger.info(colorscale)
	# getting the data specific to selected species
	base_values_ = [fig_base_values[i] for i in indices]
	nucleotide_bases_ = [fig_nucleotide_bases[i] for i in indices]
	nucleotide_bases_ = nucleotide_bases_
	nucleotide_bases2 = [item for sublist in nucleotide_bases_ for item in sublist]
	hover_values_ = [fig_hover_values[i] for i in indices]
	seq_names_ = [sp_names[i] for i in indices]
	logger.info(seq_names_)

	y = [[i]*seq_len for i in seq_names_]
	y = [item for sublist in y for item in sublist]


	sp__ = [get_key(x) for x in seq_names_]
	legend_col_ =[]
	for l in legend_col:
		tf= l[0]
		for s in sp__:
			try:
				if s in motif_species[tf].lower():
					legend_col_.append(l)
					break
			except:
				if s in motif_species[tf.split(' : ')[0]].lower(): # if multiple TF has same name and same motif
					legend_col_.append(l)
					break
	#temp = [[x[0], x[1]] for x in legend_col_]
	#temp_single = [item for sublist in temp for item in sublist]
	#colorscale_ = [[0, 'whitesmoke']] + [x for x in colorscale if x[1] in temp_single]

	#logger.info(motif_species)
	logger.info('updated legend_col_')
	logger.info(legend_col_)
	#logger.info('colorscale_')
	#logger.info(colorscale_)
	
	#logger.info('y\n{}', y)
	#print('-----------')
	#logger.info('base_values_seleceted species:\n{}',base_values_)	
	#print('-----------')
	#logger.info('colorscale inside get_figure {}',colorscale)
	#print('-----------')
	#logger.info('hover_values_seleceted species:\n{}',hover_values_[::-1])
	#logger.info('hover_values_seleceted species:\n{}',nucleotide_bases_)
	#logger.info('hover_values_seleceted species:\n{}',nucleotide_bases2)

	trace = dict(type='heatmap', z=base_values_, colorscale = colorscale, 
			 showscale=False, text=hover_values_, hoverinfo='text',zmin=0, zmax=1
			)
	data=[trace] # generating color based on the values
	data.append({'type': 'scattergl',
					'mode': 'text',
					'x': list(range(seq_len))*len(nucleotide_bases_),
					'y': y,
					 'hoverinfo':'none',
					'text': nucleotide_bases2,
					'textfont': {
						'family':'Courier New',
						'size': 15,
						
					}})
	#print(data)
	steps = [{'args': ['xaxis', {'range': [-0.5 + e, 40.5 + e]}], 'method': 'relayout', 'label':e} for e in range(seq_len-30)]
	
	sliders = [dict(active = 0,steps = steps)]
	layout = dict(sliders=sliders)
	layout['xaxis'] = {'range': [-0.5, 40.5]}
	layout['yaxis'] = {'tickfont':{'size':11}}
	layout['autosize'] = True
	#layout['width'] = 1000
	#layout['height'] = 300
	fig = dict(data=data, layout=layout)
	
	########################################## FIGURE BAR ##################################
	figbar = plotly.subplots.make_subplots(rows=3, cols=2,vertical_spacing=0.05,
                                 column_widths=[.975,.025],
                                 row_heights=[.74,.11,.15],
                                specs=[[None, {}],
                               [{"colspan": 2}, {}],
                                [{'colspan':2}, None]])
	
	
	bar = go.Heatmap(z=base_values_, colorscale=colorscale,text=hover_values_,
									 hoverinfo ='text',
									showscale=False,
									zmin=0,
									zmax=1
									 )
	
	figbar.append_trace(bar, 2, 1)
	figbar.append_trace(go.Scatter(x=[None],y=[None],
                             showlegend=False
                            ), 3,1)

	if len(legend_col_) > 0:
		for values in legend_col_:
			sm=go.Scatter(x=[None], y=[None], mode='markers',
					   marker=dict(size=10, color=values[1]),
					   legendgroup='motifs', showlegend=True, name=values[0])
			figbar.append_trace(sm, 1, 2)
		figbar.update_xaxes(showticklabels=False, row=1, col=2,showgrid=False, showline=False)
		figbar.update_yaxes(showticklabels=False, row=1, col=2, showgrid=False, showline=False)

	
	

	figbar.update_yaxes(showticklabels=False, row=2, col=1)
	figbar.update_xaxes(showticklabels=False, row=3, col=1)
	figbar.update_yaxes(showticklabels=False, row=3, col=1)

	figbar.update_layout(
		height=410,
		plot_bgcolor='white'
	)



	
	return fig, figbar, tree_diagram(seq_names_, motif_species, sp__)








