import xml.etree.ElementTree as et
from cdxml_elements import *

class CDXML_Editor(object):
	def __init__(self, cdxml_files_path=None, output_path=None):
		self.container = BOX()
		self.cdxml_path = cdxml_files_path
		self.output_path = output_path

		if self.output_path:
			self.tree = et.ElementTree()
			self.cdxml = et.Element('CDXML', attrib={'color':'3','bgcolor':'2'})
			self.tree._setroot(self.cdxml)
			self.add_colortable()
			self.page = et.SubElement(self.cdxml,'page', attrib={
				'HeightPages': '1',
				'WidthPages': '1',
				'DrawingSpace':'poster'})

	def add_colortable(self):
		self.colortable = et.SubElement(self.cdxml,'colortable')
		bg_color = et.SubElement(self.colortable,'color',attrib={'r':'1','g':'0.980','b':'0.941'}) #2 - white
		fg_color = et.SubElement(self.colortable,'color',attrib={'r':'0.200','g':'0.200','b':'0.200'}) #3 - black
		promiscuous = et.SubElement(self.colortable,'color',attrib={'r':'0.500','g':'0.500','b':'0.500'}) #4 - gray
		intermediates = et.SubElement(self.colortable,'color',attrib={'r':'0.200','g':'0.200','b':'0.800'}) #5 - blue
		target = et.SubElement(self.colortable,'color',attrib={'r':'0.800','g':'0.200','b':'0.200'}) #6 - red
		self.color_index = 6
		
	def append(self, box, arrange=""):
		self.container.append(box, arrange=arrange)
	
	def parse_cdxml(self, compound):
		path = self.cdxml_path + compound + '.cdxml'
		output = ''
		try:
			output = et.parse(path).find('.//*')
		except:
			# Extract cpd_name from path if no cdxml file exists
			output = path.split('/')[-1][:-6]
		return output	

	def get_cpd_cdxml(self, compound, id_offset=0, color="5"):
		try:
			cpd,label = self.parse_cdxml(compound)
			output = COMPOUND(cpd, label, id_offset=id_offset, color=color)
			last_id = int(output.last_id)
			return output, last_id
		except:
			cpd_name = self.parse_cdxml(compound)
			output = BOX()
			output.append(TEXT(cpd_name, color=color))
			return output, id_offset


	def add_reactants(self, reactants, previous_reactions, last_id):
		if len(previous_reactions) == 0: # First reaction
			for r_index in range(len(reactants)):
				rc, last_id = self.get_cpd_cdxml(reactants[r_index], id_offset=last_id, color="4")
				self.append(rc, arrange="right")
				if r_index < len(reactants)-1:
					self.append(TEXT('+', color="4"), arrange="right")
		else:
			for pr in previous_reactions:
				main_rs = list(set(pr.products).intersection(set(reactants)))
				box = BOX()
				for mr_index in range(len(main_rs)):
					rc, last_id = self.get_cpd_cdxml(main_rs[mr_index], id_offset=last_id, color="5")
					box.append(rc, arrange="right")
					if mr_index < len(main_rs)-1:
						box.append(TEXT('+',color="5"), arrange="right")
				pr.append(box, arrange="right")
			if len(previous_reactions) == 1:
				self.append(previous_reactions[0].container, arrange="right")
			else:
				box = BOX()
				for pr in previous_reactions:
					box.append(pr.container, arrange="below", align="right")
				box.y_reference = box.height/2
				self.append(box, arrange="right")

		return last_id


	def add_product(self, product, last_id):
		target, last_id = self.get_cpd_cdxml(product, id_offset=last_id, color="6")
		self.append(target, arrange="right")


	def add_transition(self, promiscuous_reactants, promiscuous_products, misc_products, 
						reaction_proteins="", reaction_solvents="", reaction_catalysts="", reaction_SPRESI_info=""):
		self.transition = TRANSITION(reactants=promiscuous_reactants, products=promiscuous_products, misc_products=misc_products,
										reaction_proteins=reaction_proteins, reaction_solvents=reaction_solvents,
										reaction_catalysts=reaction_catalysts, reaction_SPRESI_info=reaction_SPRESI_info)
		self.container.append(self.transition, arrange="right")

	def set_products(self, cdxml_products):
		self.products = cdxml_products

	def set_FBA(self, color_index):
		self.transition.set_color(str(color_index))

	def add_color(self, fba_value):
		if self.output_path:
			et.SubElement(self.colortable, 'color',attrib={'r':str(fba_value),'g':'0','b':'0'})
			self.color_index += 1

	def generate_file(self):	
		self.container.set_x(self.container.get_x() + 50)
		self.container.set_y(self.container.get_y() + 100)

		for element in self.container.root.findall('./'):
			self.page.append(element)
		
		self.page.set('BoundingBox', '0 0 %f %f' % (self.container.width+150, self.container.height+150))
		self.page.set('Width', str(self.container.width + 150))
		self.page.set('Height', str(self.container.height + 150))
		self.tree.write(self.output_path, encoding='UTF-8')