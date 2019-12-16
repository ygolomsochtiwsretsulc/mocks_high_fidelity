"""
Creates plots after an ETC run

A few general plots.

A plot per template.

"""


# imports
import os, sys
import numpy as n
import astropy.io.fits as fits

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

# pathes for S8 on johan's laptop
root_dir = '/home/comparat/data/4most/S8/'
template_dir = '/home/comparat/data/4most/templates/'
path_2_template = lambda tpl_name : os.path.join(template_dir, tpl_name)
path_2_catalog = os.path.join( root_dir, 'S8_Cosmology_May2nd_correct_output_ETC13_06_19_21-40-19.fits' )
output_folder = os.path.join( root_dir, 'etc_plots')
path_2_archive = os.path.join( root_dir, 'etc_plots_S8.tar.gz')


# reads the output catalog
hd = fits.open(path_2_catalog)

template = hd[1].data['TEMPLATE']
ruleset  = hd[1].data['RULESET'] 

uniq_templates_i = n.unique( template )
uniq_rulesets_i  = n.unique( ruleset  )

uniq_templates_j, uniq_rulesets_j = n.meshgrid(uniq_templates_i, uniq_rulesets_i) 

uniq_combination = n.transpose( [n.hstack((uniq_templates_j)), n.hstack((uniq_rulesets_j))]  )

def plot_tpl(uniq_c):
	print(uniq_c)
	template_name, ruleset_name =uniq_c
	data_template = fits.open(path_2_template(template_name))

	x_template = data_template[1].data['LAMBDA']
	y_template = data_template[1].data['FLUX_DENSITY']
	ok =(y_template>1e-18)

	selection = (template == template_name) & (ruleset == ruleset_name)
	N_used = len(selection.nonzero()[0])
	if N_used>1:
		print('N in catalog', N_used)
		MAG    = hd[1].data['MAG'   ][selection]
		TEXP_B = hd[1].data['TEXP_B'][selection]
		TEXP_G = hd[1].data['TEXP_G'][selection]
		TEXP_D = hd[1].data['TEXP_D'][selection]

		# one plot per template 
		# A4 figure
		fig = p.figure(0, (8.2, 11.7), frameon=False )

		# template
		fig.add_subplot(411, 
						title='template='+template_name , 
						xlabel='wavelength [Angstrom]', 
						ylabel=r'Flux [$f_\lambda$ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]',
						ylim=(( n.median(y_template[ok])/10, n.median(y_template[ok])*10 )),
						xlim=(( 3600, 9400 )),
						yscale='log' )

		p.grid()
		p.plot(x_template[ok], y_template[ok], color='black', lw=1)
		ys =  [n.min(y_template[ok]), n.max(y_template[ok])]
		p.fill_betweenx(ys, x1=[4200.0,4200.0],x2=[5000.0, 5000.0], alpha=0.1, color='b')
		p.fill_betweenx(ys, x1=[5500.0,5500.0],x2=[6700.0, 6700.0], alpha=0.1, color='g')
		p.fill_betweenx(ys, x1=[7200.0,7200.0],x2=[9000.0, 9000.0], alpha=0.1, color='r')

		# exposure time tracks
		fig.add_subplot(412,
						title='ruleset='+ ruleset_name,
						xlabel='MAG', 
						ylabel='TEXP [minutes]',
						yscale='log')

		p.plot(MAG, TEXP_B, color='r', marker='x', ls='', label=r'bright N$_{fh}$='+str(n.round(n.sum(TEXP_B)/60.,1)), rasterized=True)
		p.plot(MAG, TEXP_G, color='g', marker='+', ls='', label=r'grey N$_{fh}$='  +str(n.round(n.sum(TEXP_G)/60.,1)), rasterized=True)
		p.plot(MAG, TEXP_D, color='b', marker='^', ls='', label=r'dark N$_{fh}$='  +str(n.round(n.sum(TEXP_D)/60.,1)), rasterized=True)
		p.grid()
		p.legend( frameon=False, loc=0)

		# statistics in catalog
		fig.add_subplot(413, 
						xlabel='MAG', 
						ylabel='Counts (N)',
						yscale='log')

		p.hist(MAG, bins=n.arange(n.min(MAG)-0.1, n.max(MAG)+0.1, 0.1 ), histtype='step', lw=2)
		p.grid()
		
		fig.add_subplot(414, 
						xlabel='MAG', 
						ylabel='N x TEXP (hours)',
						yscale='log')

		p.hist(MAG, weights=TEXP_D/60., bins=n.arange(n.min(MAG)-0.1, n.max(MAG)+0.1, 0.1 ), histtype='step',label=r'bright', lw=2)
		p.hist(MAG, weights=TEXP_G/60., bins=n.arange(n.min(MAG)-0.1, n.max(MAG)+0.1, 0.1 ), histtype='step',label=r'grey' , lw=2)
		p.hist(MAG, weights=TEXP_B/60., bins=n.arange(n.min(MAG)-0.1, n.max(MAG)+0.1, 0.1 ), histtype='step',label=r'dark' , lw=2)
		p.legend( frameon=False, loc=0)

		p.grid()
		p.savefig(os.path.join(output_folder, ruleset_name+'_'+template_name[:-5]+'.png'))
		p.tight_layout()
		p.clf()
	else:
		print('not used')

for uniq_c in uniq_combination:
	plot_tpl(uniq_c)
	
os.system('tar -czf '+path_2_archive+' '+output_folder)
