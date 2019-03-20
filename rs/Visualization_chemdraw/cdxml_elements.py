import xml.etree.ElementTree as et

ARROW_LENGTH = 200
DEFAULT_Y_REFERENCE = 200
SCALE_FACTOR = 0.75

class ELEMENT(object):
	def __init__(self, el_id="-1", x_0=0, y_0=0, color="0"):
		self.el_id = el_id
		self.type = "element"
		self.x_0 = x_0
		self.y_0 = y_0
		self.color = color
		self.root = et.Element('placeholder_root')

	def get_id(self):
		return self.el_id

	def get_x(self):
		return self.x_0

	def get_y(self):
		return self.y_0

	def get_color(self):
		return self.color

	def set_id(self, new_id):
		self.el_id = new_id
		self.root.set('id', self.el_id)

	def update_position(self):
		self.root.set('p', '%f %f' % (self.x_0, self.y_0))

	def set_x(self, new_x):
		self.x_0 = new_x
		self.update_position()

	def set_y(self, new_y):
		self.y_0 = new_y
		self.update_position()	

	def set_color(self, new_color):
		self.color = new_color
		self.root.set('color', new_color)

	def append(self, element):
		try:
			self.root.append(element.root)
		except:
			self.root.append(element)

	def dump(self):
		return et.dump(self.root)

	def iter(self, query):
		return self.root.iter(query)

class ARC(ELEMENT):
	def __init__(self, el_id="-1", x_0=0, y_0=0, color="3", size=ARROW_LENGTH*SCALE_FACTOR/8, head=False):
		super(ARC, self).__init__(el_id,x_0,y_0,color)
		self.size = size
		self.head = head
		self.root = et.Element('graphic', attrib={
			'id': self.el_id,
			'GraphicType':'Arc',
			'Head3D':'%f %f 0' % (self.x_0, self.y_0),
			'Tail3D':'%f %f 0' % (self.x_0, self.y_0),
			'AngularSize':'-90',
			'HeadSize':'1000',
			'color':color
			})
		if self.head:
			self.root.set('Tail3D','%f %f 0' % (self.x_0+self.size, self.y_0))
			self.root.set('ArrowType','FullHead')
		else:
			self.root.set('Head3D','%f %f 0' % (self.x_0, self.y_0-self.size))

	def update_position(self):
		self.root.set('Head3D', '%f %f 0' % (self.x_0+self.size, self.y_0+self.size))
		self.root.set('Tail3D', '%f %f 0' % (self.x_0, self.y_0))

class ARROW(ELEMENT):
	def __init__(self, el_id="-1", x_0=0, y_0=0, color="3", width=ARROW_LENGTH*SCALE_FACTOR):
		super(ARROW, self).__init__(el_id,x_0,y_0,color)
		self.height = 0
		self.width = width
		self.y_reference = self.y_0
		self.root = et.Element('graphic', attrib={
			'BoundingBox': '%f %f %f %f' % (self.x_0+self.width, self.y_0, self.x_0, self.y_0),
			'id': self.el_id,
			'GraphicType':'Line',
			'ArrowType':'FullHead',
			'HeadSize':'1000',
			'color':color
			})

	def update_position(self):
		self.root.set('BoundingBox', '%f %f %f %f' % (self.x_0+self.width, self.y_0, self.x_0, self.y_0))

class TEXT(ELEMENT):
	def __init__(self, label, el_id="-1",x_0=0, y_0=0, color="0"):
		super(TEXT,self).__init__(el_id=el_id,x_0=x_0, y_0=y_0, color=color)
		self.label = label
		self.width = 6 * len(self.label)
		self.height = 0
		self.root = et.Element('t',attrib={
			'p':'%f %f' % (self.x_0, self.y_0),
			'Justification':'Center',
			'InterpretChemically':'no'})
		self.s = et.SubElement(self.root,'s', attrib={
			'color':color,
			'face': '96',
			'font': '3',
			'size': '10'})
		self.s.text = self.label

	def set_color(self, new_color):
		self.color = new_color
		self.s.set('color', self.color)

	def get_label(self):
		return self.label

	def set_label(self, new_label):
		self.label = new_label
		self.s.text = self.label



class BOX(ELEMENT):
	def __init__(self, width=0, height=0):
		super(BOX,self).__init__()
		self.type = "box"
		self.width = width
		self.height = height
		self.color = "3"
		self.x_reference = 0
		self.y_reference = 0
		self.isEmpty = True
		self.root = et.Element('placeholder_box')

	def append(self, element, arrange="", align="middle"):
		''' Append elements/boxes to this box '''
		if self.isEmpty: # First element in the box
			self.width = element.width
			self.height = element.height
			if hasattr(element,'type') and element.type == 'box': # Compound
				self.x_reference = element.x_reference
				self.y_reference = element.y_reference
			elif hasattr(element,'type') and element.type == 'element' and hasattr(element,'label'): # Text
				self.x_reference = element.width/2
				self.y_reference = element.y_0
				element.set_x(self.x_reference)
		elif hasattr(element, 'type') and element.type == 'element' and hasattr(element, 'label'):
			element.y_reference = 0
			element.x_reference = 0
			element.height = 0


		if arrange == 'above' or arrange == 'below' and not self.isEmpty:
			plus_symbol = TEXT('+',color="5")
			plus_symbol.set_x(max(self.x_reference, element.x_reference))
			if arrange == 'above':
				plus_symbol.set_y(element.height + 10)
			else:
				plus_symbol.set_y(-10)

			
			if align == "middle" and self.width > element.width:
				element.set_x(self.x_0 + self.width/2)
			elif align == "middle": 
				self.set_x(element.x_0 + element.width/2)
			elif align == "right":
				self.set_x(max(self.x_reference,element.x_reference) - self.x_reference)
				element.set_x(max(self.x_reference, element.x_reference) - element.x_reference)
				if arrange == 'above':
					self.y_reference = element.height
				else:
					self.y_reference = self.height

			element.append(plus_symbol)

			if arrange == 'above':
				self.set_y(element.height + 10)
			else:
				element.set_y(self.height + 10)

			self.height = self.height + element.height + 10
			self.width = max(self.width, element.width)

		elif arrange == 'right' and not self.isEmpty:
			self.x_reference = self.width + 10 + element.x_reference
			element.set_x(self.width + 10)
			element.set_y(self.y_reference - element.y_reference)
			self.width = self.width + element.width + 10
			self.height = max(self.height, element.height)

		if hasattr(element, 'type') and element.type == 'box':
			self.append_box(element)
		else:
			super(BOX,self).append(element)
		self.isEmpty = False
		
		
	def append_box(self, box_element):
		subelements = box_element.root.findall('./*')
		for obj in subelements:
			self.root.append(obj)

	def set_color(self, new_color):
		self.color = new_color
		for obj in self.root.iter():
			if obj.get('color') or obj.get('font') or obj.get('B') or obj.tag=='graphic':
				obj.set('color', self.color)

	def set_x(self, new_x):
		dx = new_x - self.get_x()
		self.x_0 = new_x
		self.update_position(dx=dx)

	def set_y(self, new_y):
		dy = new_y - self.get_y()		
		self.y_reference = self.y_reference + (new_y - self.y_0)
		self.y_0 = new_y
		self.update_position(dy=dy)

	def update_position(self, dx=0, dy=0):
		for obj in self.root.iter():
			if obj.get('p'):
				x,y = [float(n) for n in obj.get('p').split(' ')]				
				obj.set('p', '%f %f' % (x+dx, y+dy))
			elif obj.get('BoundingBox'):
				x1,y1,x2,y2 = [float(n) for n in obj.get('BoundingBox').split(' ')]
				obj.set('BoundingBox', '%f %f %f %f' % (x1+dx, y1+dy, x2+dx, y2+dy))
			elif obj.get('Head3D') and obj.get('Tail3D'):
				x,y,z = [float(n) for n in obj.get('Head3D').split(' ')]
				obj.set('Head3D', '%f %f 0' % (x+dx, y+dy))
				x,y,z = [float(n) for n in obj.get('Tail3D').split(' ')]
				obj.set('Tail3D', '%f %f 0' % (x+dx, y+dy)) 


# ARROW WITH PROMISCUOUS COMPOUNDS
class TRANSITION(BOX):
	def __init__(self, reactants, products, misc_products,
					reaction_proteins="", reaction_solvents="", reaction_catalysts="", reaction_SPRESI_info="",	color="3"):		
		super(TRANSITION,self).__init__(width=ARROW_LENGTH * SCALE_FACTOR)		
		self.color=color
		self.arrow = ARROW(color=self.color)
		self.append(self.arrow)

		self.arcs = False
		self.topArc = ARC(
			color=self.color,
			x_0=3*self.arrow.width/8,
			)
		self.bottomArc = ARC(
			color=self.color,
			x_0=5*self.arrow.width/8,
			y_0=self.arrow.width/8,
			head=True
			)		
		
		if reactants or products or misc_products:
			self.height=ARROW_LENGTH*SCALE_FACTOR/2			
			self.arcs = True
			if reactants:
				self.append(self.topArc)
			if products or misc_products:
				self.append(self.bottomArc)
		
		if reactants:
			self.reactants = TEXT(
				' + '.join(reactants), 
				x_0 = self.width/4,
				y_0 = -self.arrow.width/4,
				color = "4")
			self.append(self.reactants, align="")
		if products:
			self.products = TEXT(
				' + '.join(products), 
				x_0 = 3*self.width/4,
				y_0 = self.arrow.width/4,
				color = "4")
			self.append(self.products, align="")
		if misc_products:
			txt = '+\n'.join(misc_products)
			y_placement = self.height/2
			if products:
				txt = '+\n' + txt
				y_placement = self.height
			
			self.misc_products = TEXT(
				txt,
				x_0 = 3*self.width/4,
				y_0 = y_placement,
				color = "5")
			self.append(self.misc_products, align="")

		if reaction_proteins:
			self.reaction = TEXT(
				reaction_proteins,
				x_0 = self.arrow.width/2,
				y_0 = 10,
				color = "3")
			self.append(self.reaction)

		chemical_rxn_info = []
		if reaction_solvents:
			chemical_rxn_info += [reaction_solvents]
		if reaction_catalysts:
			chemical_rxn_info += [reaction_catalysts]
		if reaction_SPRESI_info:
			chemical_rxn_info += [reaction_SPRESI_info]
		if len(chemical_rxn_info) > 0:
			self.rxn_info = TEXT(
				'\n'.join(chemical_rxn_info),
				x_0 = self.arrow.width/2,
				y_0 = 100,
				color=self.color)
			self.append(self.rxn_info)
			self.height += 250

	def set_color(self, new_color):
		self.arrow.set_color(new_color)
		if self.arcs:
			self.topArc.set_color(new_color)
			self.bottomArc.set_color(new_color)
		if hasattr(self,'rxn_info'):
			self.rxn_info.set_color(new_color)

class COMPOUND(BOX):
	def __init__(self, cpd_xml, label_xml, id_offset=0, color="0"):
		super(COMPOUND,self).__init__()
		self.append(cpd_xml)
		self.append(label_xml, shift_down=True)
		self.cpd_id = id_offset + 1
		self.last_id = self.offset_ids(id_offset)
		self.color = color
		self.set_color(self.color)
		self.scale()
		self.set_dimensions()

	def append(self, xml_element, shift_down=False):
		if shift_down:
			for obj in xml_element.iter():
				if obj.get('p'):
					x,y = [float(n) for n in obj.get('p').split(' ')]
					obj.set('p','%f %f' % (x, y+20))
		self.root.append(xml_element)

	def scale(self):
		for obj in self.root.iter():
			if obj.get('p'):
				x,y = [float(n) for n in obj.get('p').split(' ')]
				obj.set('p','%f %f' % (x*SCALE_FACTOR, y*SCALE_FACTOR))
			if obj.get('BoundingBox'):
				x1,y1,x2,y2 = [float(n) for n in obj.get('BoundingBox').split(' ')]
				obj.set('BoundingBox','%f %f %f %f' % (x1*SCALE_FACTOR, y1*SCALE_FACTOR, x2*SCALE_FACTOR, y2*SCALE_FACTOR))

	def offset_ids(self, id_offset):		
		last_id = id_offset
		for obj in self.root.iter():
			if obj.get('id') and obj.get('id') != "-1":
				last_id = int(obj.get('id')) + int(id_offset)
				obj.set('id', str(last_id))
			if obj.get('B'):
				obj.set('B', str(int(obj.get('B')) + int(id_offset)))
				obj.set('E', str(int(obj.get('E')) + int(id_offset)))
		return str(last_id)

	def set_dimensions(self):
		max_x = 0
		max_y = 0
		x = 0
		y = 0
		for obj in self.root.iter():			
			if obj.get('p'):
				x,y = [float(n) for n in obj.get('p').split(' ')]				
			if obj.get('BoundingBox'):
				x,y,a,b = [float(n) for n in obj.get('BoundingBox').split(' ')]
			max_x = max(x, max_x)
			max_y = max(y, max_y)

		self.width = int(max_x + 30)
		self.height = int(max_y + 30)
		self.x_reference = self.width/2
		self.y_reference = self.height/2