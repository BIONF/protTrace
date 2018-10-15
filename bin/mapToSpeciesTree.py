import os, sys
import dendropy
import configure
import colourizeTree

### Module to read in the decay parameters file
### Calculate traceability for every species in the species tree
### Colourize the species tree based on traceability

def main(protId, config_file):
	rootDir = os.getcwd()

	prot_id = protId
	prot_config = configure.setParams(config_file)
	species_id = prot_config.species

	work_dir = prot_config.path_work_dir + '/' + prot_id
	os.chdir(work_dir)

	cache_dir = prot_config.path_cache

	nexus_file = 'nexus_' + prot_id + '.nexus'

	# Save the species name into a variable taxonset
	trees = dendropy.TreeList.get_from_path(prot_config.species_tree, "newick")
	taxonset = []
	for element in trees.taxon_namespace:
		taxonset.append(str(element).replace("'", "").replace(" ", "_"))

	generateNexusFile(nexus_file,taxonset,prot_config)
	colourizeTree.main(nexus_file, prot_config.hamstr_oma_tree_map, prot_id, prot_config.species_tree, prot_config.plot_figtree, prot_config.species_MaxLikMatrix, species_id, cache_dir,prot_config.fas_score)

	os.chdir(rootDir)

def generateNexusFile(nexus_file,taxonset,prot_config):
	print '##### Generating nexus file #####'
	fnew = open(nexus_file, 'w')
	fnew.write('#NEXUS\nbegin taxa;\n\tdimensions ntax=%s;\n\ttaxlabels' %len(taxonset))
	for taxa in taxonset:
		fnew.write('\n\t' + taxa.replace(' ', '_'))
	fnew.write('\n;\nend;\n\nbegin trees;\n\ttree tree_1 =\n')
	f = open(prot_config.species_tree).read()
	fnew.write(f)
	fnew.write('\nend;\n\nbegin figtree;\n')
	fnew.write('\tset appearance.backgroundColorAttribute="User Selection";\n\tset appearance.backgroundColour=#-1;\n\tset appearance.branchColorAttribute="User Selection";\n\tset appearance.branchLineWidth=3.0;\n\tset appearance.foregroundColour=#-16777216;\n\tset appearance.selectionColour=#-2144520576;\n\tset branchLabels.colorAttribute="User Selection";\n\tset branchLabels.displayAttribute="bootstrap";\n\tset branchLabels.fontName="Times New Roman";\n\tset branchLabels.fontSize=28;\n\tset branchLabels.fontStyle=1;\n\tset branchLabels.isShown=true;\n\tset branchLabels.significantDigits=4;\n\tset layout.expansion=0;\n\tset layout.layoutType="POLAR";\n\tset layout.zoom=1100;\n\tset nodeBars.barWidth=4.0;\n\tset nodeLabels.colorAttribute="User Selection";\n\tset nodeLabels.displayAttribute="Node ages";\n\tset nodeLabels.fontName="sansserif";\n\tset nodeLabels.fontSize=14;\n\tset nodeLabels.fontStyle=0;\n\tset nodeLabels.isShown=false;\n\tset nodeLabels.significantDigits=4;\n\tset polarLayout.alignTipLabels=false;\n\tset polarLayout.angularRange=0;\n\tset polarLayout.rootAngle=0;\n\tset polarLayout.rootLength=100;\n\tset polarLayout.showRoot=false;\n\tset radialLayout.spread=0.0;\n\tset rectilinearLayout.alignTipLabels=false;\n\tset rectilinearLayout.curvature=0;\n\tset rectilinearLayout.rootLength=100;\n\tset scale.offsetAge=0.0;\n\tset scale.rootAge=1.0;\n\tset scale.scaleFactor=1.0;\n\tset scale.scaleRoot=false;\n\tset scaleAxis.automaticScale=true;\n\tset scaleAxis.fontSize=8.0;\n\tset scaleAxis.isShown=false;\n\tset scaleAxis.lineWidth=1.0;\n\tset scaleAxis.majorTicks=1.0;\n\tset scaleAxis.origin=0.0;\n\tset scaleAxis.reverseAxis=false;\n\tset scaleAxis.showGrid=true;\n\tset scaleAxis.significantDigits=4;\n\tset scaleBar.automaticScale=true;\n\tset scaleBar.fontSize=10.0;\n\tset scaleBar.isShown=true;\n\tset scaleBar.lineWidth=1.0;\n\tset scaleBar.scaleRange=0.0;\n\tset scaleBar.significantDigits=4;\n\tset tipLabels.colorAttribute="User Selection";\n\tset tipLabels.displayAttribute="Names";\n\tset tipLabels.fontName="Times New Roman";\n\tset tipLabels.fontSize=18;\n\tset tipLabels.fontStyle=1;\n\tset tipLabels.isShown=true;\n\tset tipLabels.significantDigits=4;\n\tset trees.order=false;\n\tset trees.orderType="increasing";\n\tset trees.rooting=false;\n\tset trees.rootingType="User Selection";\n\tset trees.transform=false;\n\tset trees.transformType="cladogram";')
	fnew.write('\nend;')

