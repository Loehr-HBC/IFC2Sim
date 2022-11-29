# Nodes for IFC2SIM
bl_info = {
    "name": "IFC2SIM",
    "author": "Felix Loehr",
    "version": (0, 9, 20220223),
    "blender": (2, 80, 0),
    "location": "Node Editor > IFC2SIM",
    "description": "Nodes to support the complete workflow from IFC4 to EnergyPlus, including analysis",
    "warning": "",
    "doc_url": "",
    "category": "Node",
}
import bpy
from bpy.types import (NodeTree, Node, NodeSocket)
from functools import reduce # we reduce the filters onto the objects in NodeObjectsFilterLive
from mathutils import Vector # needed in the createBoundContains_test factory
from textwrap import TextWrapper # wraps our error/info text.
import readeso, os # NOTE: 'os' is already imported by 'bpy.path' as '_os'

#from .props import (PG_Objects)

# NOTE: we use BEGIN and END for wrapping/collapsing code in the KATE-texteditor

# END ##########################################################################
###                            !!! CONSTANTS !!!                             ###
#######################################80############################### BEGIN #
#   some things are treated differently between versions. (e.g. missing icons)
_B_VERSION, _B_VERSIONs = bpy.app.version, bpy.app.version_string
_B_VERS =   _B_VERSIONs[:3]
#   more readable version of the keys shortened by readeso.short_names
_dct = {'sol_diffuse'       :'Diffuse Solar',               ### site keys
        'sol_beam'          :'Direct Solar',
        'dry_bulb'          :'Dry-Bulb Temperature',
        'rel_hum'           :'Relative Humidity (Outside)',
        'w_dir'             :'Wind Direction',
        'w_speed'           :'Wind Speed',
        'air_rel_hum'       :'Relative Humidity (Zone)',    ### zone keys
        'temp_air'          :'Air Temperature',
        'temp_op'           :'Operative Temperature',
        'q_heat_equip'      :'Equipment Heating',
        'q_heat_sens'       :'Sensible Heating',
        'q_heat_ideal_sens' :'Sensible Heating (ideal)',
        'q_heat_recovery'   :'Sensible Heat Recovery',
        'p_win_sol_trans'   :'Insolation',
        'inf_ach'           :'Infiltration (Air-Change per hour)'}
# END ##########################################################################
###                            !!! FUNCTIONS !!!                             ###
#######################################80############################### BEGIN #
### Link zones that were exported, by starting from their nodes
#   NOTE: compares only matching polygons and with fixed precisions
def link_exported_zones():
    ### fetch the "EnVi-Network" node-tree
    # NOTE: the VI-Suite itself only supports one tree, so we stick to that
    for tree in bpy.data.node_groups:
        if tree.bl_idname=="EnViN": break
    else: return
    ### fetch the nodes
    all_zones = {}
    for node in tree.nodes:
        if node.bl_idname == "No_En_Net_Zone":  # "Zone" nodes
            # we may just setdefault as we may filter afterwards.
            all_zones.setdefault(node.zone, []).append(node)
    ### gather polygons for comparision
    all_polygons = dict()
    for zone in list(all_zones.keys()):# explicit, so we can kick bad items
        if zone not in bpy.data.objects:    del all_zones[zone]; continue   # broken
        node = all_zones[zone][0]# just grab any node. (guaranteed at least one)
        skts = (skt.name for skt in node.outputs if skt.bl_idname=="So_En_Net_Bound")
        if  not skts:                       del all_zones[zone]; continue   # no interzone boundaries
        # gather information
        all_polygons[zone] = dict()
        obj = bpy.data.objects[zone]
        co_verts = [v.co for v in obj.data.vertices]
        for skt in skts:
            poly  = obj.data.polygons[int(skt.split("_")[-2])]
            verts = np.array(sorted((co_verts[v] for v in poly.vertices)))
            area  = poly.area
            mid   = poly.center
            all_polygons[zone].setdefault(
                # corners       0.1 m**2   0.1m center
                f"{len(verts)} {area:.1f} {mid.to_tuple(1)}",
                list() # payload -- must be list as vectors are unhashable
                ).append((skt, [area, mid, verts])) # key-value-tuple
        # now sort and cast into OrderedDict to speed up comparision later
        for k in all_polygons[zone]: # OrderedDict[skt]=[area,mid,verts,skt]
            all_polygons[zone][k] = OrderedDict(sorted(all_polygons[zone][k],
                key=lambda amvs:(amvs[1][0],*amvs[1][1])))
    ### Compare polygons for matching vertices (if all matches, its a match)
    #   compare 'new' to 'visited' ones. that means we need'nt pop from dict.
    visited_polygon_zones, all_matches = set(), set()
    for zone_a in all_polygons:
        poly_info_a,     del_b2 = all_polygons[zone_a], set()
        for zone_b in visited_polygon_zones:
            poly_info_b, del_a2 = all_polygons[zone_b], set()
            for key_p in poly_info_a: # we have a key acting as bin, so ...
                if key_p not in poly_info_b: continue   # ... skip ahead
                matches, del_a = set(), set()
                for skt_a in poly_info_a[key_p]:
                    (    area_a, mid_a, verts_a) = poly_info_a[key_p][skt_a]
                    for skt_b in poly_info_b[key_p]:
                        (area_b, mid_b, verts_b) = poly_info_b[key_p][skt_b]
                        # we need check area (we're probably better already):
                        if min(area_a, area_b)/max(area_a, area_b) < 0.9: break
                        # ... then the centers
                        if (mid_a-mid_b).magnitude > 0.01: break
                        # still here? take a look at the vertices
                        if np.abs(verts_a - verts_b).max() < 0.005:
                            # found a match ... should we check for a better
                            # match or is the first one allways ok? Should
                            # we test, not based on absolute 5mm dist,
                            # but based on dist/area?
                            matches.add((zone_a, skt_a, skt_b, zone_b))
                            # we can remove the entry since we no longer iter
                            del poly_info_b[key_p][skt_b]
                            del_a.add(skt_a)
                            break # got our match, so skip ahead
                # cleanup
                for skt_a in del_a: del poly_info_a[key_p][skt_a]
                all_matches.update(matches)
                if not poly_info_a[key_p]:  del_a2.add(key_p)
                if not poly_info_b[key_p]:  del poly_info_b[key_p]
            for key_p in del_a2:            del poly_info_a[key_p]
            if not poly_info_b:             del_b2.add(zone_b)
        ## finally remove depleted zones from stack
        for zone_b in del_b2: visited_polygon_zones.remove(zone_b)
        if poly_info_a: # if we're blank, there's nothing to append
            visited_polygon_zones.add(zone_a)# dont forget to append ...
    ### apply the matches that we found.
    # NOTE: Linking is EXTREMELY slow (400s vs. 0.02s for all the rest).
    #       This appears to be caused by excessive iteration over node.links
    #       within the node.update-function. node.links is a getter over
    #       node_tree.links and link-changes fire every node.update.
    for (zone_a, skt_a, skt_b, zone_b) in all_matches:
        nodes_a, nodes_b = all_zones[zone_a], all_zones[zone_b]
        if len(nodes_a)+len(nodes_b)>2:
            nodes_a.sort(key=lambda n:n[0].location)
            nodes_b.sort(key=lambda n:n[0].location)
            print("Multiple zone-nodes using the same zone. We'll use the"
                +f" left-most ones:\n ('{nodes_a[0].name}' and '{nodes_b[0].name}')"
                +f"Please remove the rest:\n {[n.name for n in nodes_a + nodes_b]}")
            all_zones[zone_a], all_zones[zone_b] = [nodes_a[0]],[nodes_b[0]]
        for node_a, node_b in zip(nodes_a, nodes_b):
            links, rechts = sorted(((node_a, skt_a), (node_b, skt_b)),
                                key=lambda n:n[0].location)
            upd_a, upd_b = node_a.update, node_b.update # try overriding the
            node_a.update = node_b.update = lambda s:None # updates for speed
            node_a.id_data.links.new(   links[0].outputs[links[ 1]],
                                        rechts[0].inputs[rechts[1]])
            node_a.update, node_b.update = upd_a, upd_b # restore updates
            break # we may only support ONE link.

### fix a bug with the volume-calculation within the VI-Suite
#   TODO: rework, so it may be constrained to certain node-trees &/ objects
def set_correct_zone_volume():        # fixes a bug with the volumes in vi-suite
    depsgraph = bpy.context.evaluated_depsgraph_get()
    zonen = set()

    for tree in bpy.data.node_groups:
        if tree.bl_idname=="EnViN":                     # "EnVi-Network" node-trees
            for node in tree.nodes:
                if node.bl_idname == "No_En_Net_Zone":  # "Zone" nodes
                    if node.zone in zonen:              # already visited zone
                        node.volcalc = '0'
                    elif node.zone in bpy.data.objects:
                        node.volcalc = '0'

                        obj = bpy.data.objects[node.zone]
                        bm = bmesh.new()
                        bm.from_object(obj, depsgraph)
                        bm.transform(obj.matrix_world)
                        bmesh.ops.triangulate(bm, faces=bm.faces[:])
                        obj['auto_volume'] = bm.calc_volume()
                        bm.free()

# END ##########################################################################
###                            !!! UTILITIES !!!                             ###
#######################################80############################### BEGIN #
### Dangerous stuff: Monkey-patching private members in private objects 2.9+ ###
################################################################################
_view_layer_update = bpy.ops._BPyOpsSubModOp._view_layer_update
_view_layer_update_dummy = lambda context:None
def makeUpdateDummy():
    bpy.ops._BPyOpsSubModOp._view_layer_update = _view_layer_update_dummy
def makeUpdateValid():
    bpy.ops._BPyOpsSubModOp._view_layer_update = _view_layer_update
def RunWithoutViewlayerUpdate(function, *args, **kwargs):
    try:# we must ensure we won't crash before restoring the function
        bpy.ops._BPyOpsSubModOp._view_layer_update = _view_layer_update_dummy
        return function(*args, **kwargs)# 'finally' happens before the 'return'
    except Exception as e: raise e      # and before the Exception is re-raised
    finally: bpy.ops._BPyOpsSubModOp._view_layer_update = _view_layer_update
#####
def _copy_settings(from_this, to_that, propnames=set()):
    for prop in propnames:
        try: setattr(to_that, getattr(from_this, prop))
        except: pass

# END ##########################################################################
###                            !!! FACTORIES !!!                             ###
#######################################80############################### BEGIN #
### BEGIN Pointer-PropertyGroups (for creating Collections holding References) #
### Creates new Pointer-PropertyGroup and registeres it.
def createPointerGroup(name, ptr_type, ptr_alias, ptr_update=None, flip=False):
    ptr, ptr_alias = ('pointer', ptr_alias)[::[1,-1][ptr_alias and flip]]
    settings = {'__annotations__':{ptr:bpy.props.PointerProperty(type=ptr_type, update=ptr_update)}}
    if ptr_alias: settings[ptr_alias] = property(
        lambda s:s.__getattribute__(ptr), lambda s,p:s.__setattr__(ptr, p))
    cls = type(name, (bpy.types.PropertyGroup,), settings)
    bpy.utils.register_class(cls); return cls
### use it                      class_name      pointer_type    alias (no-update)
PG_Objects = createPointerGroup("PG_Objects", bpy.types.Object,  "object")
PG_Collections = createPointerGroup("PG_Collections", bpy.types.Collection,  "collection")
#PG_Meshes  = createPointerGroup("PG_Meshes",  bpy.types.Mesh,    "mesh"  )
#PG_Scenes  = createPointerGroup("PG_Scenes",  bpy.types.Scene,   "scene" )
################################################################################
def createPointerGroup_getter_setter(item_type):
    """Creates a factory for defining getters and setters for PointerGroups"""
    def _makeGetterSetter(col_name):
        get_items = lambda self:[ptr.pointer for ptr in getattr(self,col_name)]
        def set_items(self, items):
            collection = getattr(self,col_name); collection.clear()
            for item in items:
                if type(item)==item_type: collection.add().pointer=item
        #del_items= lambda self:getattr(self,col_name).clear()# not yet used
        return get_items, set_items
    return _makeGetterSetter # factory for factories of getters/setters
PG_Objects_getter_setter = createPointerGroup_getter_setter(bpy.types.Object)
### END Pointer-PropertyGroups (for creating Collections holding References) #
# we might want simple collections to store ints, etc, but frankly
# to just store names, bpy.types.PropertyGroup is enough
def createByName_getter_setter(source, col_name):
    get_items = lambda self:[source.get(p.name) for p in getattr(self,col_name)]
    def set_items(self, items):
        collection = getattr(self,col_name); collection.clear()
        for item in items:# assumes hasattr(item, "name") to be true
            if item in source: collection.add().name = item.name
    #del_items= lambda self:getattr(self,col_name).clear()# not yet used
    return get_items, set_items
# this only works for AABBs
#def createBoundContains_test(obj):# creates a test for bounding boxes
    #A, G = obj.bound_box[0], obj.bound_box[6]# get diagonal corners
    #A, G = obj.matrix_world @ Vector(A), obj.matrix_world @ Vector(B)# to world
    #a, g = zip(*[(min(z), max(z)) for z in zip(A,G)])# sort min and max values
    #return lambda vec:all((a[i]<=vec[i]<=g[i] for i in [1,2,3]))
# however if i transform the tested point ... Still, we wouldn't get crossings.
#def createBoundContains_test(obj):# creates a test for bounding boxes
    #A, G = obj.bound_box[0], obj.bound_box[6]# get diagonal corners
    #A, G = obj.matrix_world @ Vector(A), obj.matrix_world @ Vector(B)# to world
    #a, g = zip(*[(min(z), max(z)) for z in zip(A,G)])# sort min and max values
    #return lambda vec:all((a[i]<=vec[i]<=g[i] for i in [1,2,3]))
################################################################################
class PG_Filter(bpy.types.PropertyGroup):# work with filters
    # NOTE: This is only the right-hand-part of a lambda function taking 'obj'.
    eval_line:bpy.props.StringProperty(name="lambda obj:", description='Lambda '
        'expression. Evaluate to True if to keep the provided object "obj".')
    invert: bpy.props.BoolProperty(description="Invert this filter")
    @property
    def filter_func(self): # quick-and-dirty. Chains errors if unfunctional.
        if not self.eval_line: return lambda obj:True
        return eval("lambda obj:"+["{}","not ({})"][self.invert].format(self.eval_line))
bpy.utils.register_class(PG_Filter)
class PG_Filter_B(bpy.types.PropertyGroup):# work with filters
    invert: bpy.props.BoolProperty(  description="Invert this filter")
    string: bpy.props.StringProperty(description="string for comparision")
    ### filters
    filter_type:bpy.props.EnumProperty(name="", description="Main filter type.",
        items=[ ("flt_name",    "Name",      "Filter by Object-Name"          ),
                ("flt_IFCType", "IFC-Type",  "Filter by IFC-Type"             ),
                ("flt_ZoneType","ZoneType",  "Filter by type of zone"         ),
                ("flt_origin",  "Origin",    "Filter by Object-Origin"        ),
                ("flt_vertex",  "Vertex",    "Filter by Vertex-Location"      ),
                ("flt_center",  "Center",    "Filter by Center-of-Geometry"   ),
                ("flt_slice",   "Slice",
                    "Filter by position in list as python slice (start,stop,step)"),
                ("eval_line",   "Expression","Filter by any python-expression"),
                ])
    _items=[("WEST",  "West",  "Largest X <"), ("EAST", "East", "Smallest X >"),
            ("SOUTH", "South", "Largest Y <"), ("NORTH","North","Smallest Y >"),
            ("BELOW", "Below", "Largest Z <"), ("ABOVE","Above","Smallest Z >"),
            ("DISTANCE", "Distance", "Distance"), ("BOX", "Box", "Box")     ]
    flt_center:bpy.props.EnumProperty( name="", description="Filter by Center.",
        items=_items)
    flt_origin:bpy.props.EnumProperty( name="", description="Filter by Vertex-Coordinate",
        items=_items)
    flt_vertex:bpy.props.EnumProperty( name="", description="Filter by Vertex-Coordinate",
        items=_items)
    flt_distance:bpy.props.EnumProperty(name="", description="Filter by Vertex-Coordinate",
        items=[ ("LT",      "Less",     "Less Than"         ),
                ("GT",      "Greater",  "Greater Than"      ),
                ("BETWEEN", "Between",  "Between A and B"   ),
                ("OUTSIDE", "Outside",  "Outside A and B"   ), ])
    flt_box     :bpy.props.EnumProperty(name="", description="Filter by Vertex-Coordinate",
        items=[ ("BETWEEN", "Between",  "Between A and B"   ),
                ("OUTSIDE", "Outside",  "Outside A and B"   ),
                ("OVERLAP", "Overlap",  "Overlap with box"  ), ])
    flt_slice: bpy.props.IntVectorProperty(size=3, default=(0,0,1), name="slice",
        description="Slice-indices: [start:stop:step] or, if start==stop [start::step]")
    vec_lsw: bpy.props.FloatVectorProperty(description="Lower-South-West corner")
    vec_tne: bpy.props.FloatVectorProperty(description="Higher-North-East corner")
    vec_to : bpy.props.FloatVectorProperty(description="Point for Distance")
    float_min:bpy.props.FloatProperty(description="Minimum float to compare to")
    float_max:bpy.props.FloatProperty(description="Maximum float to compare to")
    bool_flp:bpy.props.BoolProperty()
    #
    flt_name  : bpy.props.EnumProperty( name="", description="Filter name by.",
        items=[ ("STARTSWITH", "Begin",         "object.name.startswith"),
                ("ENDSWITH",   "End",           "object.name.endswith"  ),
                ("CONTAINS",   "Contains",      "string in object.name" ),
                ("CONTAINED",  "Contained By",  "object.name in string" ),  ])
    flt_IFCType:bpy.props.EnumProperty(name="", description="Filter by Type.",
        items=[ ("IfcSpace",    "IfcSpace", "IfcSpace; general (no subtype)" ),
                ("IfcWindow",   "IfcWindow","IfcWindow; general (no subtype)"),
                ("IfcDoor",     "IfcDoor",  "IfcDoor;   general (no subtype)"),
                # openings are searched two-step
                ("Opening",       "Opening",       "IfcOpeningElement; general"),
                ("Opening/Window","Opening/Window","IfcOpeningElement; Window"),
                ("Opening/Door",  "Opening/Door",  "IfcOpeningElement; Door"),
                ])
    flt_ZoneType:bpy.props.EnumProperty(name="",description="Filter by ZoneType.",
        items=[ ("BATH",    "Bath",         "(string in 'LongName')"),
                ("HALL",    "Hall",         "(string in 'LongName')"),
                ("LIVING",  "Living-room",  "(string in 'LongName')"),
                ("SLEEPING","Sleeping-room","(string in 'LongName')"),
                ])
    # NOTE: This is only the right-hand-part of a lambda function taking 'obj'.
    eval_line:bpy.props.StringProperty(name="lambda obj:", description='Lambda '
        'expression. Evaluate to True if to keep the provided object "obj".')
    @property
    def filter_func(self): # quick-and-dirty. Chains errors if unfunctional.
        if not self.eval_line: return lambda obj:True
        return eval("lambda obj:"+["{}","not ({})"][self.invert].format(self.eval_line))
    def filtered(self, objects):
        if   self.filter_type == "flt_slice"    :# try fetching by python slice
            s,e,w = self.flt_slice
            if s==e: e=None
            if w==0: w=1
            return list(objects)[s:e:w]#[slice(*self.flt_slice)]
        if   self.filter_type == "flt_ZoneType" :# try fetching by zone type
            # to gather the relevant possible substrings
            _matches_ = {"BATH":{"BATH","BAD","WC","KLO","TOILET","TOILETTE","RESTROOM"},
                         "HALL":{"HALL","FLUR","ENTRANCE","EINGANG","GANG"},
                       "LIVING":{"LIVING","WOHN"},
                     "SLEEPING":{"BED","SLEEPING","SCHLAF"},
                         }.get(self.flt_ZoneType, {self.flt_ZoneType})
            def flt(obj):# True if either LongName or Name contains a match
                from blenderbim.bim.ifc import IfcStore
                from ifcopenshell.api.attribute.data import Data
                oprops = obj.BIMObjectProperties
                if oprops.ifc_definition_id not in Data.products:
                    Data.load(IfcStore.get_file(), oprops.ifc_definition_id)
                d = Data.products[oprops.ifc_definition_id]
                keys = sorted([(ln["name"]=="LongName", ln["value"].upper())
                                for ln in d if ln["name"]in{"LongName","Name"}],
                        key=lambda ln:ln[0])
                return any( map(keys[0][1].__contains__, _matches_)) or \
                    any(    map(keys[1][1].__contains__, _matches_))
            return filter(flt, objects)
        elif self.filter_type == "flt_IFCType"  :# try fetching by ifc-type
            if self.flt_IFCType.startswith("Ifc"):
                flts = [lambda obj:obj.name.startswith(self.flt_IFCType+"/")]
            elif self.flt_IFCType.startswith("Opening"):# opening-elements
                flts = [lambda obj:obj.name.startswith("IfcOpeningElement")]
                # the following handles special openings like windows or doors
                subtypes = {"Opening/Window" : ["window", "fenster"],
                            "Opening/Door"   : ["door","tuer","tür"],}
                if self.flt_IFCType in subtypes:
                    flts.append(lambda obj:any(map(obj.name.lower().__contains__,
                        subtypes[self.flt_IFCType])))
            return reduce(lambda x,y:filter(y,x), flts, objects)
        elif self.filter_type == "flt_name"     :# try fetching by name
            if self.flt_name=="STARTSWITH":flt=lambda s:s.startswith(self.string)
            if self.flt_name=="ENDSWITH"  :flt=lambda s:s.endswith(  self.string)
            if self.flt_name=="CONTAINS"  :flt=lambda s:self.string in s
            if self.flt_name=="CONTAINED" :flt=lambda s:s in self.string
            return filter(lambda obj:flt(obj.name), objects)
        elif self.filter_type == "eval_line"    :# try fetching by any expressio
            return filter(self.filter_func, objects)
        ### we got the other filters already. Gather geometric data.
        if   self.filter_type == "flt_origin"   :# try fetching by origin
            _obj = lambda obj:[obj.location]
            _typ = self.flt_origin
        elif self.filter_type == "flt_center"   :# try fetching by center
            def _obj(obj):
                if obj.type == "MESH": return [obj.location]
                verts = obj.data.vertices
                return [sum((v.co for v in verts[1:]),verts[0].co)/len(verts)]
            _typ = self.flt_center
        elif self.filter_type == "flt_vertex"   :# try fetching by vertex
            if  self.flt_vertex in {"WEST",  "EAST" }:  key=lambda c:c.x
            if  self.flt_vertex in {"SOUTH", "NORTH"}:  key=lambda c:c.y
            if  self.flt_vertex in {"BELOW", "ABOVE"}:  key=lambda c:c.z
            def _obj(obj):
                if obj.type != "MESH": return[obj.location]
                return sorted((v.co for v in obj.data.vertices),key=key)
#            flt = lambda obj:_obj()
            _typ= self.flt_vertex
            if self.flt_vertex in {"DISTANCE", "BOX"}:
                def _obj(obj):
                    if obj.type != "MESH": return [obj.location]
                    return [v.co for v in obj.data.vertices]
        if _typ=="DISTANCE":
            _dist = lambda obj:[(v - self.vec_to).magnitude for v in _obj(obj)]
            if self.flt_distance=="LT":
                flt=lambda obj: self.float_min> min(_dist(obj))
            if self.flt_distance=="GT":
                flt=lambda obj: self.float_max< max(_dist(obj))
            if self.flt_distance=="BETWEEN":
                flt=lambda obj:(self.float_min<=min(_dist(obj)) and
                                self.float_max>=max(_dist(obj)) )
            if self.flt_distance=="OUTSIDE":
                flt=lambda obj:(self.float_min> min(_dist(obj)) or
                                self.float_max< max(_dist(obj)) )
        elif _typ=="BOX":
            _in=lambda vec:(self.vec_to.x<=vec.x<=self.vec_to.x and
                            self.vec_to.y<=vec.y<=self.vec_to.y and
                            self.vec_to.z<=vec.z<=self.vec_to.z )
            if self.flt_box=="BETWEEN":
                flt=lambda obj:all(    _in(v)for v in _obj(obj))
            if self.flt_box=="OUTSIDE":
                flt=lambda obj:all(not _in(v)for v in _obj(obj))
            if self.flt_box=="OVERLAP":# only for vertices
                #flt=lambda obj:any(_in(v)for v in obj)# OVERLAP or BETWEEN
                def flt(obj): # OVERLAP only
                    vecs = _obj(obj)
                    ret  = _in(vecs[0])
                    for v in vecs[1:] :
                        if ret!=_in(v): return True
                    return False
        # All options for 'vectors'
        elif _typ=="WEST":  flt=lambda obj:_obj(obj)[-1].x<self.float_min
        elif _typ=="EAST":  flt=lambda obj:_obj(obj)[ 0].x>self.float_max
        elif _typ=="SOUTH": flt=lambda obj:_obj(obj)[-1].y<self.float_min
        elif _typ=="NORTH": flt=lambda obj:_obj(obj)[ 0].y>self.float_max
        elif _typ=="BELOW": flt=lambda obj:_obj(obj)[-1].z<self.float_min
        elif _typ=="ABOVE": flt=lambda obj:_obj(obj)[ 0].z>self.float_max
        return filter(self.filter_func, objects)
    def draw(self, layout):
        row = layout.row(align=True)
        row.prop(self, "filter_type")   # get type of filter (name, location, etc)
        row.prop(self, self.filter_type)# subtype (min, min-x, etc)
        if self.filter_type=='flt_name':
            row = layout.row(align=True)
            row.prop(self, "invert")    # invert result
            row.prop(self, "string")    # string to compare to
        #if self.filter_type=='flt_IFCType':
        #if self.filter_type=='flt_ZoneType':
        if self.filter_type in {'flt_origin', 'flt_center', 'flt_vertex'}:
            if   getattr(self,self.filter_type)=='DISTANCE':
                row = layout.row(align=True)# vector to compare to
                row.prop(self, "vec_to")
                row = layout.row(align=True)# float(s) to compare to
                if self.flt_distance!="GT": row.prop(self, "float_min")
                if self.flt_distance!="LT": row.prop(self, "float_max")
            elif getattr(self,self.filter_type)=='BOX':
                # vector to compare to
                layout.row().prop(self, "vec_lsw", text="LSW")
                layout.row().prop(self, "vec_tne", text="TNE")
            elif getattr(self,self.filter_type)in{'EAST','NORTH','ABOVE'}:
                row = layout.row(align=True)# float to compare to
                #row.label(text=" of:")
                row.prop(self, "float_max")
            elif getattr(self,self.filter_type)in{'WEST','SOUTH','BELOW'}:
                row = layout.row(align=True)# float to compare to
                #row.label(text=" of:")
                row.prop(self, "float_min")
        #
bpy.utils.register_class(PG_Filter_B)
################################################################################
# which settings to use for creating or updating HVAC, occupancy, etc.
# setting to NONE removes said Technode. 'GUESS' was meant to setup by educated
# guess, but honestly, we can just default our Dummy-Nodes to our educated guess
# and push that for update. => None|Custom(Update)
class PG_Technode(bpy.types.PropertyGroup):
    #--_none_custom_guess = [('NONE',"None","dont create"),
    #--    ('CUSTOM',"Custom","Custom Setup, with the values you provide"),
    #--    ('GUESS', "Guess", "Take an educated guess based on the type of the"
    #--     " first zone. Then set up accordingly.")]
    # just store, if to use a node at all. If so, on first update we'll init the
    technodes : bpy.props.BoolVectorProperty(size=8) # input-node with our guess
    # the second four bits are to store the previous state, as to know changes
    ###### draw utility
    def draw(self, node, layout, index):
        layout.context_pointer_set('node', node)
        for i,txt in enumerate(("HVAC","Infiltration","Equipment","Occupancy")):
            if not i%2: row = layout.row(align=True)
            row2 = row.row(align=True)
            row2.prop(self, "technodes", index=i, toggle=True, text=txt)
            row2.alert=self.technodes[i]!=self.technodes[i+4]
            op = row2.operator("ifc2sim.update_technodes", text="", icon="FILE_REFRESH")
            op.index, op.mode = index, txt
bpy.utils.register_class(PG_Technode)
# ---------------------------
# | Zones_name | update all |
# ---------------------------
# | HVAC|()    | Infiltr|() |
# | Occupan|() | Equipme|() |
# ---------------------------

################################################################################


# END ##########################################################################
###                             !!! UI-LISTS !!!                             ###
#######################################80############################### BEGIN #
#universal ui-list class for pointer-collections
class POINTERS_UL_general(bpy.types.UIList):
    # The draw_item function is called for each item of the collection that is visible in the list.
    #   data is the RNA object containing the collection,
    #   item is the current drawn item of the collection,
    #   icon is the "computed" icon for the item (as an integer, because some objects like materials or textures
    #   have custom icons ID, which are not available as enum items).
    #   active_data is the RNA object containing the active property for the collection (i.e. integer pointing to the
    #   active item of the collection).
    #   active_propname is the name of the active property (use 'getattr(active_data, active_propname)').
    #   index is index of the current item in the collection.
    #   flt_flag is the result of the filtering process for this item.
    #   Note: as index and flt_flag are optional arguments, you do not have to use/declare them here if you don't
    #         need them.
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        # draw_item must handle the three layout types... Usually 'DEFAULT' and 'COMPACT' can share the same code.
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            is_active = index==getattr(active_data,active_propname)
            if hasattr(item, "pointer") and hasattr(item.pointer, "name"):
                layout.prop(item.pointer, "name", text="", emboss=False,
                icon=["RIGHTARROW","TRIA_RIGHT"][is_active and self.layout_type=='DEFAULT'])
                            #icon="TRIA_RIGHT", text="", emboss=False)
            elif item.name:
                layout.prop(item, "name", text="", emboss=False, icon_value=icon)
            else:
                layout.label(text="", translate=False, icon_value=icon)
        # If we guarantee 'pointer', showing it would be suitable for 'COMPACT'.
        #elif self.layout_type in {'COMPACT'}: layout.prop(item, "pointer")
        # 'GRID' layout type should be as compact as possible (typically a single icon!).
        # there is however no sensible way of doing that for pointers
        #elif self.layout_type in {'GRID'}: #layout.label(text='', icon_value=icon)

#simple ui-list class for filter-collections
class FILTERS_UL_simple(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        # draw_item must handle self.layout_type in {'DEFAULT','COMPACT','GRID'}
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            is_active = index==getattr(active_data, active_propname)
            layout.prop(item, "eval_line", text="", emboss=False,
                icon=["RIGHTARROW","TRIA_RIGHT"][is_active and self.layout_type=='DEFAULT'])
            #if item.name:
                #layout.prop(item, "name", text="", emboss=False, icon_value=icon)
        elif self.layout_type=="GRID": layout.label(text=str(index))

#simple ui-list class for technode-collections
class TECHNODE_UL_simple(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            sockets= [inp for inp in data.inputs if inp.bl_idname=="SktObjects"]
            is_active = index==getattr(active_data, active_propname)
            if index<len(sockets):
                col = layout.column(align=True)
                col.context_pointer_set('node', data)
                # name and update-all
                row = col.row(align=True)
                row.prop(sockets[index], "name", text='', emboss=False,
                    icon=["RIGHTARROW","DOWNARROW_HLT"][is_active and self.layout_type=='DEFAULT'])
                row.alert=item.technodes[:4]!=item.technodes[4:8]
                op = row.operator("ifc2sim.update_technodes", text="", icon="FILE_REFRESH")
                op.index, op.mode = index, "ALL"
                return # wird jetzt unten gezeichnet
                # Technodes (only draw for the active item)
                if is_active: item.draw(data, col, index)
            else:
                layout.alert=True
                layout.label("Missmatch with Sockets")


# END ##########################################################################
###                             !!! SOCKETS !!!                              ###
#######################################80############################### BEGIN #
# NodeTree for IFC2SIM nodes
class IFC_SIM_Tree(NodeTree):
    '''A custom node tree for translating IFCs to vi-suite simulations.'''
    bl_idname = 'IFC2SIM_Tree'          # Optional identifier string.
    bl_label = 'IFC-to-Simulation Tree' # Label for nice name display
    bl_icon = 'NODETREE'                # Icon identifier

    # Update is called e.g. after a link was added.
    # => Use that to remove ill-defined links? Just do that in node.insert_link (which is called earlier)
    #def update(self): pass

# Socket for ...
class MyCustomSocket(NodeSocket):
    '''Custom node socket type'''
    bl_idname = 'CustomSocketType'
    bl_label = "Custom Node Socket"

    # Enum items list
    my_items = (
        ('DOWN', "Down", "Where your feet are"),
        ('UP', "Up", "Where your head should be"),
        ('LEFT', "Left", "Not right"),
        ('RIGHT', "Right", "Not left"),
    )

    my_enum_prop: bpy.props.EnumProperty(
        name="Direction",
        description="Just an example",
        items=my_items,
        default='UP',
    )

    # Optional function for drawing the socket input value
    def draw(self, context, layout, node, text):
        if self.is_output or self.is_linked:
            layout.label(text=text)
        else:
            layout.prop(self, "my_enum_prop", text=text)

    # Socket color
    def draw_color(self, context, node):
        return (1.0, 0.4, 0.216, 0.5)


# Socket for Object-Pointer-Collections
# QUASI FERTIG.
#   Soll noch löschen-nach-objekt-picker eingebaut werden?
#   Soll eine Python-Property die Objekte zurück-geben?
#       (test=>kann die dann dargestellt werden?)
class SktObjects(NodeSocket):
    '''Socket type for multiple Objects'''
    bl_idname = 'SktObjects'
    bl_label = "Objects Socket"

    ### UPDATE-FUNCTIONS (for the properties below)
    # Move reference 'object' to 'objects_collection'. If already contained,
    # the 'index' is updated to its position and the 'object' is NOT cleared.
    def upd_ptr_add(self, context):
        if self.ptr_add: self.ptr_add=False;return  # Reset button
        if not self.object: return                  # Nothing to do.
        len_oc = len(self.objects_collection)
        # we need to loop over it anyways. By looping backwards we can clean up
        for idx in range(len_oc): # empty elements (appear when objects get deleted)
            ptr = self.objects_collection[len_oc-1-idx]
            if   ptr.object==self.object: self.index = len_oc-1-idx; return
            elif ptr.object==None: self.objects_collection.remove(len_oc-1-idx)
        ptr = self.objects_collection.add()
        ptr.object, ptr.name, self.object = self.object, self.object.name, None
        self.index = len(self.objects_collection)-1
    # Pop reference from 'objects_collection' into 'object'. This way we allow
    def upd_ptr_del(self, context): # for quick re-adding in case of mistakes.
        if self.ptr_del: self.ptr_del=False;return  # Reset button
        if not self.objects_collection: return      # Nothing to do
        len_oc = len(self.objects_collection)
        loc = self.index if 0<=self.index<len_oc else len_oc-1
        self.object, self.index = self.objects_collection[loc].object, loc-1
        self.objects_collection.remove(loc)
    # ensure 'to_index' follows (thereby support moving up/down by one directly)
    def upd_idx(self, context):
        if self.to_index!=self.index: self.to_index=self.index
    # move reference within 'objects_collection' up/down/to-index
    def upd_mv(self, context):
        if self.index==self.to_index: return        # Nothing to do
        loc, len_oc = self.to_index, len(self.objects_collection)
        if   loc < 0: loc = len_oc-1 # Wrap around from top to bottom
        elif loc > len_oc-1: loc = 0 # Wrap around from bottom to top
        self.objects_collection.move(self.index, loc)   # move
        self.index = self.to_index = loc            # update index
    ### PROPERTIES
    object: bpy.props.PointerProperty(type=bpy.types.Object)
    objects_collection: bpy.props.CollectionProperty(
        name="Objects", type=PG_Objects)
    index: bpy.props.IntProperty(min=-1, default=-1, update=upd_idx)
    # function-buttons (cheaply fake non-undo-able operators)
    ptr_add: bpy.props.BoolProperty(update=upd_ptr_add)# links object in to list
    ptr_del: bpy.props.BoolProperty(update=upd_ptr_del)# remove object from list
    to_index: bpy.props.IntProperty(min=-1, default=-1, update=upd_mv)# move it.
    force_draw:bpy.props.IntProperty(min=-1,max=1,default=0)# cheat the draw-call

    # dereferenced objects
    _objects=property(lambda s:[ptr.pointer for ptr in s.objects_collection])
    @property
    def objects(self):
        if self.is_output and hasattr(self.node, "objects"): # check in the node
            return self.node.objects    # This way, we may daisy-chain getters.
        elif self.is_output or not self.is_linked: return self._objects# just us
        else: return self.links[0].from_socket.objects# fetch from linked output
    @objects.setter
    def objects(self, objs):
        self.objects_collection.clear()
        for obj in objs:# We only add objects, as other values would fail.
            if type(obj)==bpy.types.Object: # Also skips blank items
                self.objects_collection.add().pointer = obj

    ### Drawing the socket. Provides an interface including an Object-picker.
    def draw(self, context, layout, node, text):
        # if force_draw is -1 close, if it is +1 open
        if max((self.is_output or self.is_linked)-self.force_draw, 0)or not self.show_expanded:
            layout.label(text=text)# + "({})".format(len(self.objects)))
            # if self.objects comes from a longer chain, it's expensive ...
            # ... even more so, if it was filtered => Sorry, but no numbers here
        else: self._draw(context, layout, node, text)
    #   let's outsource the main segment to a helper, so nodes may use it
    def _draw(self, context, layout, node, text):
        layout=layout.column()
        row=layout.row(align=True)
        row.prop(self, 'object',  text='')
        row.prop(self, 'ptr_add', text='', icon="PLUS")
        row.prop(self, 'ptr_del', text='', icon="X")# pseudo-random ID
        layout.template_list("POINTERS_UL_general", "random"+text+node.name,
            self, "objects_collection", self, "index", rows=2, maxrows=5)
        if self.objects_collection: # display index and support moving items
            layout.prop(self, 'to_index', text="move to")# (=> update moves)

    # Socket color. Needed function. Display the color 'NodeSocketObject' uses.
    def draw_color(self, context, node):# we cannot access it directly so guess.
        return (1, 0.65, 0.4, 1)#self.is_linked*0.5+0.5)# note: result is darker

# Socket for Object-Pointer-Collections
# QUASI FERTIG.
#   Soll noch löschen-nach-objekt-picker eingebaut werden?
#   Soll eine Python-Property die Objekte zurück-geben?
#       (test=>kann die dann dargestellt werden?)
#   NodeSocketCollection ist WEISS
class SktCollections(NodeSocket):
    '''Socket type for multiple Objects'''
    bl_idname = 'SktCollections'
    bl_label = "Collections Socket"

    ### UPDATE-FUNCTIONS (for the properties below)
    # Move reference 'collection' to 'collections_collection'. If already contained,
    # the 'index' is updated to its position and the 'collection' is NOT cleared.
    def upd_ptr_add(self, context):
        if self.ptr_add: self.ptr_add=False;return  # Reset button
        if not self.collection: return              # Nothing to do.
        len_oc = len(self.collections_collection)
        # we need to loop over it anyways. By looping backwards we can clean up
        for idx in range(len_oc): # empty elements (appear when objects get deleted)
            ptr = self.collections_collection[len_oc-1-idx]
            if   ptr.collection==self.collection: self.index = len_oc-1-idx; return
            elif ptr.collection==None: self.collections_collection.remove(len_oc-1-idx)
        ptr = self.collections_collection.add()
        ptr.pointer, ptr.name = self.collection, self.collection.name
        self.collection = None
        self.index = len(self.collections_collection)-1
    # Pop reference from 'collections_collection' into 'collection'. This way we allow
    def upd_ptr_del(self, context): # for quick re-adding in case of mistakes.
        if self.ptr_del: self.ptr_del=False;return  # Reset button
        if not self.collections_collection: return  # Nothing to do
        len_oc = len(self.collections_collection)
        loc = self.index if 0<=self.index<len_oc else len_oc-1
        self.collection = self.collections_collection[loc].pointer
        self.collections_collection.remove(loc)
        self.index = loc-1
    # ensure 'to_index' follows (thereby support moving up/down by one directly)
    def upd_idx(self, context):
        if self.to_index!=self.index: self.to_index=self.index
    # move reference within 'collections_collection' up/down/to-index
    def upd_mv(self, context):
        if self.index==self.to_index: return        # Nothing to do
        loc, len_oc = self.to_index, len(self.collections_collection)
        if   loc < 0: loc = len_oc-1 # Wrap around from top to bottom
        elif loc > len_oc-1: loc = 0 # Wrap around from bottom to top
        self.collections_collection.move(self.index, loc)   # move
        self.index = self.to_index = loc            # update index
    ### PROPERTIES
    collection: bpy.props.PointerProperty(type=bpy.types.Collection)
    collections_collection: bpy.props.CollectionProperty(type=PG_Collections)
    index: bpy.props.IntProperty(min=-1, default=-1, update=upd_idx)
    icon_nr: bpy.props.IntProperty(min=-1, max=8, default=-1,
        description="Color for the icon, if any. max=8, 0=white, -1=no icon")
    # function-buttons (cheaply fake non-undo-able operators)
    ptr_add: bpy.props.BoolProperty(update=upd_ptr_add)# links object in to list
    ptr_del: bpy.props.BoolProperty(update=upd_ptr_del)# remove object from list
    to_index: bpy.props.IntProperty(min=-1, default=-1, update=upd_mv)# move it.
    force_draw:bpy.props.IntProperty(min=-1,max=1,default=0)# cheat the draw-call

    # dereferenced collections
    _collections=property(lambda s:[ptr.pointer for ptr in s.collections_collection])
    @property
    def collections(self):
        if self.is_output and hasattr(self.node, "collections"): # check in the node
            return self.node.collections # This way, we may daisy-chain getters.
        elif self.is_output or not self.is_linked: return self._collections# just us
        else: return self.links[0].from_socket.collections# fetch from linked output
    @collections.setter
    def collections(self, colls):
        self.collections_collection.clear()
        for col in colls:# We only add objects, as other values would fail.
            if type(col)==bpy.types.Collection: # Also skips blank items
                self.collections_collection.add().pointer = col

    ### Drawing the socket. Provides an interface including an Object-picker.
    def draw(self, context, layout, node, text):
        icon="MATERIAL" if _B_VERS<"2.9" else [ # 28x has no OUTLINER_COLLECTION -icon
            "COLLECTION_COLOR_0{}".format(self.icon_nr),
            "OUTLINER_COLLECTION", None][ self.icon_nr<=0 + self.icon_nr==-1]
        # if force_draw is -1 close, if it is +1 open
        if max((self.is_output or self.is_linked)-self.force_draw, 0):
            layout.label(text=text)# + "({})".format(len(self.collections)))
            # if self.collections comes from a longer chain, it's expensive ...
            # ... even more so, if it was filtered => Sorry, but no numbers here
        else: self._draw(context, layout, node, text)
    #   let's outsource the main segment to a helper, so nodes may use it
    def _draw(self, context, layout, node, text):
        layout=layout.column()
        row=layout.row(align=True)
        row.prop(self, 'collection',text='')
        row.prop(self, 'ptr_add',   text='', icon="PLUS")
        row.prop(self, 'ptr_del',   text='', icon="X")# pseudo-random ID
        layout.template_list("POINTERS_UL_general", "random"+text+node.name,
            self, "collections_collection", self, "index", rows=2, maxrows=5)
        if self.collections_collection: # display index and support moving items
            layout.prop(self, 'to_index', text="move to")# (=> update moves)

    # Socket color. Needed function. Display the color 'NodeSocketObject' uses.
    def draw_color(self, context, node):# we cannot access it directly so guess.
        #return (1, 0.4, 0.65, 1)#self.is_linked*0.5+0.5)# note: result is darker
        return (1, 1, 1, 1)#self.is_linked*0.5+0.5)# note: result is darker


# END ##########################################################################
###                              !!! NODES !!!                               ###
#######################################80############################### BEGIN #
# Mix-in class for all custom nodes in this tree type.
# Defines a poll function to enable instantiation.
class IFC_SIM_Tree_Node:
    @classmethod
    def poll(cls, ntree):
        return ntree.bl_idname == IFC_SIM_Tree.bl_idname

# standard socket types
sockettypes = [ # These sockets are predefined
    ('NodeSocketBool',             "bool"      ),
    ('NodeSocketColor',            "color"     ),

    ('NodeSocketFloat',            "Float"     ),
    *[('NodeSocketFloat'+s, s)  for s in ("Angle", "Distance", "Factor",
                                        "Percentage", "Time", "Unsigned")],

    ('NodeSocketGeometry',         "Geometry"  ),
    ('NodeSocketImage',            "Image"     ),

    *[('NodeSocket'+s, s)       for s in ("Int", "IntFactor",
                                        "IntPercentage", "IntUnsigned")],

    ('NodeSocketObject',           "Object"    ),
    ('NodeSocketShader',           "Shader"    ),

    # ('NodeSocketStandard',         "Standard"  ), # => undefined types
    ('NodeSocketString',           "String"    ),

    ('NodeSocketVector',           "Vector"    ),
    *[('NodeSocketVector'+s, s) for s in ("Acceleration",
        "Direction", "Euler", "Translation", "Velocity", "XYZ")],

    ('NodeSocketVirtual',          "Virtual"   ) ]

# Derived from the Node base type.
class DummyNode(Node, IFC_SIM_Tree_Node):
    # === Basics ===
    '''A dummy node. Not registered. It's just for reference''' # Description
    bl_idname = 'DummyNode-type' # Identifier. If not given, python class name.
    bl_label = "Dummy Node" # Label. Displayed if no function nor explicit label
    bl_icon = 'SOUND' # Icon identifier
    #bl_description = "" # ? Description?
    #bl_height_min, ..._max, ..._default # hoehe des Nodes (min, max, default)
    #bl_width_min,  ..._max, ..._default # weite des Nodes (min, max, default)

    #dimensions (readonly. groesse des Nodes)
    #height, width, width_hidden # hoehe, breite, breite wenn minimiert

    # === Custom Properties ===
    # These work just like custom properties in ID data blocks
    # Extensive information can be found under
    # http://wiki.blender.org/index.php/Doc:2.6/Manual/Extensions/Python/Properties
    my_string_prop: bpy.props.StringProperty()
    my_float_prop: bpy.props.FloatProperty(default=3.1415926)


# Derived from the Node base type.
class SocketTesterNode(Node, IFC_SIM_Tree_Node):
    # === Basics ===
    '''A dummy node which carries all sockets to just test them''' # Description
    bl_idname = 'SocketTesterNode' # Identifier. If not given, python class name.
    bl_label = "Custom Node" # Label for nice name display
    bl_icon = 'SOUND' # Icon identifier

    # === Optional Functions ===
    # Initialization function, called when a new node is created.
    # This is the most common place to create the sockets for a node, as shown below.
    # NOTE: this is not the same as the standard __init__ function in Python, which is
    #       a purely internal Python method and unknown to the node system!
    def init(self, context):
        for k, n in sockettypes: # testing all the standard socket types
            self.inputs.new( k, n+"_in" )
            self.outputs.new(k, n+"_out")
        # BUGFIX: the direction is spawned with (0,0,0), which blocks the widget
        self.inputs["Direction_in"].default_value = [1,0,0] # this un-sticks it

    # Update on node-graph-changes (nodes added/removed or links added/removed)
    def update(self): pass

    def socket_value_update(self, context): # update on socket changes
        print(context==bpy.context, context)

    #def insert_link(link): # what is that?
    # Node-Link-collection for muting nodes: self.internal_links

    # Copy function to initialize a copied node from an existing one.
    def copy(self, node):
        print("Copying from node ", node)

    # Free function to clean up on removal.
    def free(self):
        print("Removing node ", self, ", Goodbye!")

    # Additional buttons displayed on the node.
    def draw_buttons(self, context, layout):
        #layout.label(text="Node settings")
        #layout.prop(self, "my_float_prop")
        layout.prop(self.inputs["Direction_in"], "default_value")
        layout.prop(self.inputs["Direction_in"], "default_value", expand=True)


    # Detail buttons in the sidebar.
    # If this function is not defined, the draw_buttons function is used instead
    def draw_buttons_ext(self, context, layout):
        layout.prop(self, "my_float_prop")
        # my_string_prop button will only be visible in the sidebar
        layout.prop(self, "my_string_prop")

    # Optional: custom label
    # Explicit user label overrides this, but here we can define a label dynamically
    def draw_label(self):
        return "I am a custom node"


class SocketValueTesterNode(Node, IFC_SIM_Tree_Node):
    '''A dummy node for testing the value and type of sockets'''
    bl_idname, bl_label, bl_icon = ['SocketValueTesterNode']*2 + ['SOUND']
    def init(self, context):
        # NodeSocketVirtual seems to be the one at bottom of the GroupInputs
        for k, n in [('NodeSocketVirtual', 'Virtual')]:
            self.inputs.new( k, n+"_in" )
            self.outputs.new(k, n+"_out")

    def draw_buttons(self, context, layout):
        show = []
        if self.inputs[0].links:
            ln = self.inputs[0].links[0]
            show.append(("Input", ln.from_node, ln.from_socket, ln))
        if self.outputs[0].links:
            ln = self.outputs[0].links[0]
            show.append(("Output", ln.to_node, ln.to_socket, ln))
        for title, node, socket, link in show:
            box = layout.box()
            box.label(text=title)
            box2 = box.box()
            box2.label(text="Node")
            box2.prop(node, "name")
            box2.prop(node, "bl_idname")
            box3 = box.box()
            box3.label(text="Socket")
            box3.prop(socket, "name")
            box3.prop(socket, "bl_idname")
            if hasattr(socket, 'default_value'):
                box3.label(text=str(type(socket.default_value)))
                box3.label(text="%r"%socket.default_value)
                box3.prop(socket, 'default_value')
        if not show:
            layout.label(text="Nothing connected")

class SocketNode(Node, IFC_SIM_Tree_Node):
    '''A dummy node for testing the value and type of sockets'''
    bl_idname, bl_label, bl_icon = ['SocketNode']*2 + ['SOUND']
    def init(self, context):
        # NodeSocketVirtual seems to be the one at bottom of the GroupInputs
        for k, n in [('SktObjects', 'SktObjects')]:
            self.inputs.new( k, n+"_in" )
            self.outputs.new(k, n+"_out")

################################################################################
#   OBJECTS : Input, Filter, Join, Display                                     #
################################################################################
# FIXME: Currently we provide filters as their own collection and we just eval
# provided python-expressions. We may pull the filtering interface into the node
# right here, however. Then each filter needn't store but a StringProperty which
# contains the filtering function. References to other filters could be supported
# by passing the collection as 'filters_collection' when evaluating. This allows
# constructing complex filters here; the add-function parses them into an item.
# In case of expressions we could add an "eval_line" here and test wether they
# are dangerous. The 'filters_collection' becomes a named string-collection.
# FIXME: We do not support intermediate access to the objects. Each filter only
# gets one object, not the "currently available" ones. => no comparing between.
class NodeObjectsFilterLive(Node, IFC_SIM_Tree_Node):
    '''Filter objects by certain properties'''
    bl_idname, bl_label, bl_icon = 'NodeObjectsFilterLive', "Objects Filter Live", 'CUBE'
    def init(self, context):
        self.inputs.new( 'SktObjects', "Objects")
        self.outputs.new('SktObjects', "remaining")

    ### UPDATE-FUNCTIONS (for the properties below)
    def upd_flt_add(self, context):# Adds a filter, but doesn't set it up
        if self.flt_add: self.flt_add=False;return  # Reset button
        self.filters_collection.add()
        self.index = len(self.filters_collection)-1
    def upd_flt_del(self, context):# removes filter: current if given else last.
        if self.flt_del: self.flt_del=False;return  # Reset button
        if not self.filters_collection: return                 # Nothing to do
        len_oc = len(self.filters_collection)
        loc = self.index if 0<=self.index<len_oc else len_oc-1
        self.filters_collection.remove(loc); self.index = loc-1
    # ensure 'to_index' follows (thereby support moving up/down by one directly)
    def upd_idx(self, context):
        if self.to_index!=self.index: self.to_index=self.index
    # move filter within 'filters_collection' up/down/to-index
    def upd_mv(self, context):
        if self.index==self.to_index: return        # Nothing to do
        loc, len_oc = self.to_index, len(self.filters_collection)
        if   loc < 0: loc = len_oc-1 # Wrap around from top to bottom
        elif loc > len_oc-1: loc = 0 # Wrap around from bottom to top
        self.filters_collection.move(self.index, loc)   # move
        self.index = self.to_index = loc            # update index

    ### PROPERTIES
    filters_collection: bpy.props.CollectionProperty(type=PG_Filter)
    index: bpy.props.IntProperty(min=-1, default=-1, update=upd_idx)# currently
    # function-buttons (cheaply fake non-undo-able operators)
    flt_add: bpy.props.BoolProperty(update=upd_flt_add)# adds filter to the list
    flt_del: bpy.props.BoolProperty(update=upd_flt_del)# remove filter from list
    to_index: bpy.props.IntProperty(min=-1, default=-1, update=upd_mv)# move it.

    # dereferenced filters
    filters=property(lambda s:[flt.filter_func for flt in s.filters_collection])
    objects=property(lambda s:list(reduce(lambda x,y:filter(y,x), s.filters,
            s.inputs[0].objects)))# filters are applied in order. => faster

    def draw_buttons(self, context, layout):
        row=layout.row(align=True)
        row.label(text="active Filters ({})".format(len(self.filters_collection)))
        row.prop(self, 'flt_add', text='', icon="PLUS")
        row.prop(self, 'flt_del', text='', icon="X")
        if self.filters_collection:# use a pseudo-random ID to not sync layouts
            layout.template_list("FILTERS_UL_simple", "random"+self.name,
                self, "filters_collection", self, "index", rows=2, maxrows=5)
            # display index and support moving items
            layout.prop(self, 'to_index', text="move to")# (=> update moves)

# FIXME: This enforces updating objects ... (should we check in each draw call?)
class NodeObjectsFilter(Node, IFC_SIM_Tree_Node):
    '''Split objects by individual filters and output them'''
    bl_label, bl_icon = "Objects Filter", 'FILTER'
    # we use this prop as identifier to ensure we dont accidentally share a
    id_string : bpy.props.StringProperty()# ui-list layout
    def init(self, context):
        self.inputs.new( 'SktObjects', "Objects")
        self.id_string = str(id(self))

    ### UPDATE-FUNCTIONS (for the properties below)
        # WARNING: WE MUST KEEP FILTER AND SOCKET IN SYNC ON OUR OWN! Outside
        # meddling may break stuff! There's no simple way of overriding any prop
        # of the item. We cant PointerProperty the node, nor the socket, nor can
        # we update props within any collection-item to refer to the socket.
        # We may not even use 'path_from_id' to gather the item (=>unsupported).
    def upd_flt_add(self, context):# Adds a filter&socket, but doesn't set it up
        if self.flt_add: self.flt_add=False;return  # Reset button
        # TODO: here we should do cleanups for "unbreaking" stuff. It might work
        # by defining a string in the item, which keeps the identifier of the
        # socket and by searching where those dont fit. The identifier is set on
        # socket-definition as its name and should be unique (maybe hash(id)?)
        self.filters_collection.add().name="undef_objects"  # filter flt
        self.outputs.new('SktObjects',"undef_objects")      # socket skt
        self.index = len(self.filters_collection)-1
    def upd_flt_del(self, context):# removes filter: current if given else last.
        if self.flt_del: self.flt_del=False;return  # Reset button
        if not self.filters_collection: return      # Nothing to do
        len_oc = len(self.filters_collection)
        loc = self.index if 0<=self.index<len_oc else len_oc-1
        self.filters_collection.remove(loc)
        self.outputs.remove(self.outputs[loc])
        self.index = loc-1
    # ensure 'to_index' follows (thereby support moving up/down by one directly)
    def upd_idx(self, context):
        if self.to_index!=self.index: self.to_index=self.index
    # move filter/socket within 'filters_collection'/'outputs' up/down/to-index
    def upd_mv(self, context):
        if self.index==self.to_index: return        # Nothing to do
        loc, len_oc = self.to_index, len(self.filters_collection)
        if   loc < 0: loc = len_oc-1 # Wrap around from top to bottom
        elif loc > len_oc-1: loc = 0 # Wrap around from bottom to top
        self.filters_collection.move(self.index, loc)   # move filter
        self.outputs.move(self.index, loc)              # move filter
        self.index = self.to_index = loc                # update index
    # update the objects of an output
    def upd_update_obj(self, context):
        if self.update_obj: self.update_obj=False;return# Reset button
        objs = self.inputs[0].objects
        flt, skt = self.filters_collection[self.index], self.outputs[self.index]
        skt.objects = flt.filtered(objs)
    def upd_update_all(self, context):
        if self.update_all: self.update_all=False;return# Reset button
        objs = self.inputs[0].objects
        for flt, skt in zip(self.filters_collection, self.outputs):
            skt.objects = flt.filtered(objs)

    ### PROPERTIES
    # add checkmark wether we need to update our results (which index?)
    filters_collection: bpy.props.CollectionProperty(type=PG_Filter_B)
    index: bpy.props.IntProperty(min=-1, default=-1, update=upd_idx)# currently
    def name_dummy_set(self, name):
        idx = self.index
        self.filters_collection[idx].name = self.outputs[idx].name = name
    name_dummy:bpy.props.StringProperty(name="Name", set=name_dummy_set,
        get=lambda s:s.outputs[s.index].name if s.outputs and s.index<=len(s.outputs)else "",
        description="shared name-property for both, socket and filter")
    # function-buttons (cheaply fake non-undo-able operators)
    flt_add: bpy.props.BoolProperty(update=upd_flt_add)# adds filter to the list
    flt_del: bpy.props.BoolProperty(update=upd_flt_del)# remove filter from list
    to_index: bpy.props.IntProperty(min=-1, default=-1, update=upd_mv)# move it.
    update_obj:bpy.props.BoolProperty(update=upd_update_obj)# single output
    update_all:bpy.props.BoolProperty(update=upd_update_all)# all outputs object

    ## dereferenced filters / we can't use that directly.
    #filters=property(lambda s:[flt.filter_func for flt in s.filters_collection])
    #objects=property(lambda s:list(reduce(lambda x,y:filter(y,x), s.filters,
    #        s.inputs[0].objects)))# filters are applied in order. => faster

    def draw_buttons(self, context, layout):
        # the show_options would hide the buttons, except for the extbuttons
        #layout.prop(self, "show_options", toggle=True)
        row=layout.row(align=True)
        row.label(text="active Filters ({})".format(len(self.filters_collection)))
        row.prop(self, 'flt_add', text='', icon="PLUS")
        row.prop(self, 'flt_del', text='', icon="X")
        if self.filters_collection:# use a pseudo-random ID to not sync layouts
            layout.template_list("POINTERS_UL_general", self.id_string, #FILTERS_UL_simple
                self, "outputs", self, "index", rows=2, maxrows=5)      #filters_collection
            # display index and support moving items
            layout.prop(self, 'to_index', text="move to")# (=> update moves)
        if 0<=self.index<=len(self.filters_collection):
            box = layout.box()
            flt = self.filters_collection[self.index]
            skt = self.outputs[self.index]
            #box.prop(self, "name_dummy")
            flt.draw(box.column(align=True))
        row=layout.row(align=True)
        if 0<=self.index<=len(self.filters_collection):
            row.prop(self, 'update_obj', text='Update',     toggle=True)
        row.prop(self, 'update_all', text='Update all', toggle=True)

# TODO/FIXME: If we don't care about order we could simply update the SktObjects
#   'objects'-property to return from all connected inputs instead of the first.
# NOTE: If we move from the pretty interface (with NodeSocketVirtual=unconnected)
#   to the simple one (last one is free) code becomes more efficient and simpler
class NodeObjectsJoin(Node, IFC_SIM_Tree_Node):
    '''Concatenate individual objects to socket'''
    bl_idname, bl_label, bl_icon = 'NodeObjectsJoin', 'Objects Join', 'CUBE'
    unique:     bpy.props.BoolProperty(description="Drop duplicates")
    keep_order: bpy.props.BoolProperty(description="Keep order", default=True)
    def init(self, context):
        self.outputs.new('SktObjects', "Objects")
        self.inputs.new( 'SktObjects', "Objects")
        # using this socket for unconnected/dummy is more pretty, but also slow
        #self.inputs.new( 'NodeSocketVirtual', "")# we would need safety
    # This is called directly after a link was added and only on those two nodes
    def insert_link(self, link):# after it returns the NodeTree.update is called
        if link.to_socket==self.inputs[-1]:
            len_in = len(self.inputs)
            a = self.inputs.new( 'SktObjects', "Objects") # add new input-socket
            a.force_draw = -1
            ### we don't use other socket-types any longer, ...
            #self.inputs.move(len_in, len_in-1)            # correct sorting
            # since we cannot replace the sockets of a link, add a new one instead
            #self.id_data.links.new(link.from_socket, self.inputs[-2])
            #self.id_data.links.remove(link)# don't forget to remove the old link
    # Update on node-graph-changes (nodes added/removed or links added/removed)
    def update(self): # removing would need to be done here?
        for skt in self.inputs:
            if not skt.is_linked and not skt==self.inputs[-1]:
                self.inputs.remove(skt)
    @property
    def objects(self): # return the [ordered] [unique] concatenation of objects
        ret = sum((skt.objects for skt in self.inputs if skt.name),[])
        if self.unique:
            if not self.keep_order: return list(set(ret))
            ret2, tmp = [], set()
            for obj in ret:
                if obj not in tmp:
                    ret2.append(obj)
                    tmp.add(obj)
            return ret2
        return ret

    def draw_buttons(self, context, layout):
        layout.prop(self, "unique", text="Unique", toggle=True)
        if self.unique:
            layout.prop(self, "keep_order", text="keep order", toggle = True)

# TODO: Add support for 'from-filter'-input
class NodeObjectsInput(Node, IFC_SIM_Tree_Node):
    '''Input objects individually or add assortments to socket.'''
    bl_idname, bl_label, bl_icon = 'NodeObjectsInput', 'Objects Input', 'CUBE'
    ### UPDATE-FUNCTIONS (for the properties below)
    # Add objects to 'outputs[0].objects_collection', depending on mode. We dont
    # respect the previous ordering, but at least we clear broken elements.
    def upd_ptr_add(self, context):
        if self.ptr_add: self.ptr_add=False;return      # Reset button
        collection = self.outputs[0].objects_collection # our collection
        # let's grab our objects
        if  self.mode=="SELECTED": objects = bpy.context.selected_objects
        # skip duplicates. frankly let's crash the order. We present it sorted:
        objects_old = {ptr.pointer for ptr in collection if ptr.pointer!=None}
        objects_old.update(objects)
        objects_old.discard(None)# discards blank entries too. Safe for sorting
        collection.clear() # clear the old objects
        for obj in sorted(objects_old, key=lambda obj:obj.name):
            ptr = collection.add()
            ptr.pointer, ptr.name = obj, obj.name
        self.outputs[0].index = len(collection)-1
    # Remove objects from 'outputs[0].objects_collection', depending on mode.
    # Respects ordering and clears broken elements.
    def upd_ptr_del(self, context):
        if self.ptr_del: self.ptr_del=False;return  # Reset button
        collection = self.outputs[0].objects_collection
        if not collection: return                   # Nothing to do
        # let's grab our objects
        if  self.mode=="SELECTED": objects = bpy.context.selected_objects
        # we need to loop backwards, as to keep our indices valid
        len_oc, index = len(collection), self.outputs[0].index
        objects = set((*objects, None))# sets are faster, None to remove broken
        for idx in range(len_oc):
            ptr = collection[len_oc-1-idx]
            if ptr.pointer in objects:
                collection.remove(len_oc-1-idx)
                index -= len_oc-1-idx <= index# move index with the items
        self.outputs[0].index = index
    ### PROPERTIES
    mode : bpy.props.EnumProperty(items=[(b.upper(),b,"") for b in ("Selected",)])
    # function-buttons (cheaply fake non-undo-able operators)
    ptr_add: bpy.props.BoolProperty(update=upd_ptr_add)# links object in to list
    ptr_del: bpy.props.BoolProperty(update=upd_ptr_del)# remove object from list
    def init(self, context): self.outputs.new('SktObjects', "Objects")
    def draw_buttons(self, context, layout):
        self.outputs[0]._draw(context, layout, self, "Objects")
        layout.prop(self, "mode")
        row = layout.row(align=True)
        row.prop(self, "ptr_add", text="Add", icon="PLUS")
        row.prop(self, "ptr_del", text="Remove", icon="X")
# NOTE: UI-Lists do not support lists and Collections do not support getters.
#   Therefore we have to help ourselves by comparing the items and writing them
#   to our own Collection, which we display in an UI-List. => Update-Button
#   The alternative would be to provide a paginated layout all by ourselves.
class NodeObjectsViewer(Node, IFC_SIM_Tree_Node):
    '''A dummy node for testing the value and type of sockets'''
    bl_idname, bl_label = 'NodeObjectsViewer', 'Objects Display'
    def upd_update_objects(self, context):
        if self.update_objects: self.update_objects=False;return
        self.objects_collection.clear()
        for obj in self.inputs[0].objects:
            self.objects_collection.add().pointer=obj
    index:  bpy.props.IntProperty(min=-1,default=-1)
    objects_collection:bpy.props.CollectionProperty(type=PG_Objects)
    update_objects:bpy.props.BoolProperty(update=upd_update_objects)
    objects=property(lambda s:[ptr.pointer for ptr in s.objects_collection])
    def init(self, context): self.inputs.new('SktObjects', "Objects")
    def draw_buttons(self, context, layout):
        skt = self.inputs[0]
        if skt.is_linked:
            if self.objects!=self.inputs[0].objects:
                row=layout.row();row.alert=True
                row.prop(self, "update_objects", toggle=True)
            layout.template_list("POINTERS_UL_general", "random"+self.name+skt.name,
                self, "objects_collection", self, "index", rows=2, maxrows=5)

# TODO: Add support for 'from-filter'-input
class NodeCollectionsAddObjects(Node, IFC_SIM_Tree_Node):
    '''Add objects or collections to (other) collections.'''
    bl_label, bl_icon = 'Add to Collections', 'CUBE'
    ### UPDATE-FUNCTIONS (for the properties below)
    # Add inputs: one SktCollections and one NodeSocketVirtual, the latter
    # will revert-from/convert-to SktCollections or SktObjects on (un-)linking.
    def upd_ptr_add(self, context):
        if self.ptr_add: self.ptr_add=False;return      # Reset button
        skt_virt = self.inputs.new('NodeSocketVirtual', "Add these items..."  )
        skt_coll = self.inputs.new('SktCollections',    "... to those Collections"   )
        skt_coll.display_shape = ['CIRCLE', 'SQUARE', 'DIAMOND'][2] # + "_DOT"
        skt_coll.force_draw=-1
        self.index = len(self.inputs)-1
    # Remove inputs: one SktCollections and one NodeSocketVirtual, which
    # might be transformed into a SktCollections or a SktObjects socket.
    def upd_ptr_del(self, context):
        if self.ptr_del: self.ptr_del=False;return              # Reset button
        idx=int(len(self.inputs)/2-1 if self.index==-1 else self.index//2)# wrap
        skt_virt, skt_coll =self.inputs[idx*2], self.inputs[idx*2+1]#get sockets
        self.inputs.remove(skt_virt);self.inputs.remove(skt_coll)#remove sockets
        self.index=idx*2-2
    # ensure 'to_index' follows (thereby support moving up/down by one directly)
    def upd_idx(self, context):
        if self.to_index!=self.index//2: self.to_index=self.index//2
    # move filter/socket within 'filters_collection'/'outputs' up/down/to-index
    def upd_mv(self, context):
        if self.to_index==self.index//2: return  # Nothing to do
        idx, loc, len_oc =self.index//2, self.to_index, len(self.inputs)//2
        if   loc < 0: loc = len_oc-1    # Wrap around from top to bottom
        elif loc > len_oc-1: loc = 0    # Wrap around from bottom to top
        up, dwn = loc<idx, loc>idx      # we move upward
        self.inputs.move(idx*2   , loc*2+dwn)   # move virtual socket
        self.inputs.move(idx*2+up, loc*2+1  )   # move collections-socket
        self.index, self.to_index  = loc*2, loc # update index
    def init(self, context): self.id_string = str(id(self))
    ### PROPERTIES
    index:  bpy.props.IntProperty(min=-1, default=-1, update=upd_idx)# currently
    id_string : bpy.props.StringProperty()# ui-list layout
    single: bpy.props.BoolProperty(description="Only update individual sockets")
    # function-buttons (cheaply fake non-undo-able operators)
    to_index:bpy.props.IntProperty(min=-1, default=-1, update=upd_mv)# currently
    ptr_add: bpy.props.BoolProperty(update=upd_ptr_add)# links object in to list
    ptr_del: bpy.props.BoolProperty(update=upd_ptr_del)# remove object from list
    def draw_buttons(self, context, layout):
        row = layout.row(align=True)
        row.prop(self, "single", text="Current")
        row.operator_menu_enum("ifc2sim.collection_items","mode")
        row = layout.row(align=True)
        row.prop(self, "ptr_add", text="Add", icon="PLUS")
        row.prop(self, "ptr_del", text="Remove", icon="X")
        if  self.inputs: # display index and support moving items
            layout.template_list("POINTERS_UL_general", self.id_string,
                self, "inputs", self, "index", rows=2, maxrows=5)
            layout.prop(self, 'to_index', text="move to")# (=> update moves)

################################################################################
#   COLLECTIONS : Input   :TODO: Modify, SortObjectsTo, Setup                  #
################################################################################
class NodeCollectionsInput(Node, IFC_SIM_Tree_Node):
    '''Input object-collections individually or add assortments to socket.'''
    bl_idname, bl_label= 'NodeCollectionsInput', 'Collections Input'
    bl_icon = 'MATERIAL'if _B_VERS<"2.9" else 'OUTLINER_COLLECTION'# erst ab 2.9
    ### UPDATE-FUNCTIONS (for the properties below)
    # Add objects to 'outputs[0].objects_collection', depending on mode. We dont
    # respect the previous ordering, but at least we clear broken elements.
    def upd_ptr_add(self, context):
        if self.ptr_add: self.ptr_add=False;return          # Reset button
        collection = self.outputs[0].collections_collection # our collection
        # let's grab our collections
        if  self.mode=="ACTIVE":
            colls = [bpy.context.view_layer.active_layer_collection.collection]
        # skip duplicates. frankly let's crash the order. We present it sorted:
        colls_old = {ptr.pointer for ptr in collection if ptr.pointer!=None}
        colls_old.update(colls)
        colls_old.discard(None)# discards blank entries too. Safe for sorting
        collection.clear() # clear the old objects
        for col in sorted(colls_old, key=lambda col:col.name):
            ptr = collection.add(); ptr.pointer = col
            #ptr.name = col.name    # we dont want names. They show blanks.
        self.outputs[0].index = len(collection)-1
    # Remove objects from 'outputs[0].objects_collection', depending on mode.
    # Respects ordering and clears broken elements.
    def upd_ptr_del(self, context):
        if self.ptr_del: self.ptr_del=False;return  # Reset button
        collection = self.outputs[0].objects_collection
        if not collection: return                   # Nothing to do
        # let's grab our collections
        if  self.mode=="ACTIVE":
            colls = [bpy.context.view_layer.active_layer_collection.collection]
        # we need to loop backwards, as to keep our indices valid
        len_oc, index = len(collection), self.outputs[0].index
        objects = set((*objects, None))# sets are faster, None to remove broken
        for idx in range(len_oc):
            ptr = collection[len_oc-1-idx]
            if ptr.pointer in objects:
                collection.remove(len_oc-1-idx)
                index -= len_oc-1-idx <= index# move index with the items
        self.outputs[0].index = index
    ### PROPERTIES
    mode : bpy.props.EnumProperty(items=[(b.upper(),b,"") for b in ("Active",)])
    # function-buttons (cheaply fake non-undo-able operators)
    ptr_add: bpy.props.BoolProperty(update=upd_ptr_add)# links object in to list
    ptr_del: bpy.props.BoolProperty(update=upd_ptr_del)# remove object from list
    def init(self, context): self.outputs.new('SktCollections', "Collections")
    def draw_buttons(self, context, layout):
        self.outputs[0]._draw(context, layout, self, "Collections")
        layout.prop(self, "mode")
        row = layout.row(align=True)
        row.prop(self, "ptr_add", text="Add", icon="PLUS")
        row.prop(self, "ptr_del", text="Remove", icon="X")

# TODO: If we just cut whenever, we might cut and cut again=> bad idea
#   => check if we need to recut => internal copy of inputs to compare against.
#   => ALSO: Do we really want to cut on the original items? maybe backup first?
class NodeCutOpenings(Node, IFC_SIM_Tree_Node):
    bl_label = "Cut Openings"
    # which icon to choose for this? EVENT_OS ? Not SYSTEM, but maybe ...
    # MOD_BOOLEAN, SELECT_INTERSECT, OBJECT_DATAMODE, OBJECT_HIDDEN, MESH_GRID
    bl_icon  = "MOD_MESHDEFORM"#"DESKTOP"# MOD_BUILD, SNAP_FACE_CENTER,
    def init(self, context):
        # NOTE: Maybe allow handing over the objects after the cut? No priority.
        #self.outputs.new('SktObjects', "Objects")
        self.inputs.new( 'SktObjects', "Zones"  ).force_draw=-1
        self.inputs.new( 'SktObjects', "Doors"  ).force_draw=-1
        self.inputs.new( 'SktObjects', "Windows").force_draw=-1
    material_door   : bpy.props.PointerProperty(type=bpy.types.Material)
    material_window : bpy.props.PointerProperty(type=bpy.types.Material)

    def draw_buttons(self, context, layout):
        #=> this needs be an operator with an undo function! Complexity!
        #layout.context_pointer_set("node",self)# already set
        op = layout.operator("IFC2SIM.cut_openings")
        layout.prop(self, 'material_door',   text="Door")
        layout.prop(self, 'material_window', text="Window")

# TODO: Design: set materials, based on assumption.
#       We should allow excluding certain materials from being replaced
#       We should support working only on selected faces or certain indices.
#####
#       Then, we should add an operator that can handle finding contacting faces
#       handle rebuilding them to apply matching materials. but that's another story
class NodeSetMaterials(Node, IFC_SIM_Tree_Node):
    bl_label = "Set Materials" # currently we only allow simple assumptions
    bl_icon  = ["FACE_MAPS","IMAGE","MATERIAL"][0] # which icon to choose for this?
    def init(self, context):
        self.inputs.new( 'SktObjects', "Zones"  ).force_draw=-1
    material_ceiling: bpy.props.PointerProperty(type=bpy.types.Material)
    material_floor  : bpy.props.PointerProperty(type=bpy.types.Material)
    material_north  : bpy.props.PointerProperty(type=bpy.types.Material)
    material_east   : bpy.props.PointerProperty(type=bpy.types.Material)
    material_south  : bpy.props.PointerProperty(type=bpy.types.Material)
    material_west   : bpy.props.PointerProperty(type=bpy.types.Material)
    def upd_set_materials(self, context):
        if self.set_materials: self.set_materials=False; return
        # gather materials
        materials = []
        for c in ("Ceiling","Floor","North","East","South","West"):
            materials.append(getattr(self, "material_"+c.lower()))#material/None
            # we know the materials must exist at this point. Use None to skip.
        # TODO: ensure objectmode somehow
        for obj in self.inputs["Zones"].objects:
            ### we might have empty objects => do skipping
            ##if obj.data==None: print(obj.name) ;continue
            # prepare materials as needed
            mat_idxs = []
            for mat in materials:
                # we use empty inputs (=> None) for skipping
                if mat==None: mat_idxs.append(None);continue
                idx = obj.data.materials.find(mat.name)
                if idx==-1:
                    idx += len(obj.data.materials)
                    obj.data.materials.append(mat)
                mat_idxs.append(idx)
            # go through faces and assign materials as fit
            for face in obj.data.polygons:
                x,y,z = face.normal
                # get index of material-meaning to set
                idx = sorted(enumerate([z,-z,y,x,-y,-x]), key=lambda i:-i[1])[0][0]
                # set the material meant for the axis, skipping empty entries
                if mat_idxs[idx]!=None:
                    face.material_index = mat_idxs[idx]# needs objectmode!
    set_materials   : bpy.props.BoolProperty(update=upd_set_materials)

    def draw_buttons(self, context, layout):
        #=> this should be an operator with an undo function! Complexity!
        #op = layout.operator("IFC2SIM.set_materials")
        layout.prop(self, 'set_materials', text="Set Materials", toggle=True)
        for c in ("Ceiling","Floor","North","East","South","West"):
            layout.prop(self, 'material_'+c.lower(), text=c)

################################################################################
#   Setup: Hardware NOTE: HVAC, etc.
# TODO: Create interface to setup zones.
#       Honestly, if someone wants to custom-setup the nodes, it might be more
#       sensible to just create said nodes in our tree, than to recreate their
#       whole interface to our node? If something changes, it will be patching
#       itself on its own. If the ID changes, that's 1-3 lines. Etc.
#       Still: adding VI-Suite-Nodes to our tree remains "probably stupid&bad"!!
# NOTE: benoetigte funktionen:
#           load_from_tree (überträgt settings vom envi_network zu uns)
#           custom_to_tree (überträgt settings von uns zum envi_network)
#           unlink_technode (entfernt den Link zum Tech-Node, nicht ihn selbst)
#           remove_technode (entfernt den Tech-Node selbst, nicht nur den Link)
# TODO: Vereinfachen: wir können einfach das Interface weiterleiten wie beim
#   EnVi Context node. => Update müsste nur noch das hinzufügen/wegnehmen machen
#   und nicht mehr das updaten => direkt in die toggles einbauen?
# TODO: operatoren/interface für die sockets unter die ui-liste zeichnen wie bei
#   split objects by filter
class NodeTechAdd(Node, IFC_SIM_Tree_Node):
    '''Add nodes for building services (HVAC, etc) to the EnVi-Network'''
    bl_idname, bl_label, bl_icon = 'NodeTechAdd', 'Building Services', 'SETTINGS'
    # pointer to EnVi-Network NOTE:  if you pointer defines a specific subtype
    #   you cant change it from GUI, only from Python and only to said subtype.
    envi_network : bpy.props.PointerProperty(type=bpy.types.NodeTree)
    #mode : bpy.props.EnumProperty(items=[("Custom","Custom","Custom"),
        #("GUESS", "Guess", "Takes an educated guess based on the type of first"
         #" zone. Then set up accordingly.")])
    #hvac :## None/Custom/Guess + Remove/unlink-from
    def init(self, context):
        self.id_string = str(id(self))# setup id
        # set up EnVi-Network pointer (If none exists, create one)
        for nt in bpy.data.node_groups:
            if nt.bl_idname == 'EnViN':
                self.envi_network = nt; break
        else:
            self.envi_network = bpy.data.node_groups.new("EnVi-Network",'EnViN')
        ## Add Input for zones to connect to shared tech.
        #for i in range(5):
            #self.inputs.new( 'SktObjects', "Zones")
            #self.technodes.add()
    def upd_ptr_add(self, context):
        if self.ptr_add: self.ptr_add=False;return              # Reset button
        self.inputs.new( 'SktObjects', "Zones").show_expanded=False
        self.technodes.add()
        #node = self.envi_network.nodes.new(type)
    def upd_ptr_del(self, context):
        if self.ptr_del: self.ptr_del=False;return              # Reset button
        idx=len(self.inputs)-1 if self.index==-1 else self.index# wrap
        inps=[i for i,skt in enumerate(self.inputs)if skt.bl_idname=="SktObjects"]
        # gather all indices from the previous to the current/target zone-socket
        idxs = list(range(inps[idx-1]+1 if idx else 0, inps[idx]+1))# current
        self.technodes.remove(idx)
        for i in idxs[::-1]: self.inputs.remove(self.inputs[i]) #remove sockets
    # ensure 'to_index' follows (thereby support moving up/down by one directly)
    def upd_idx(self, context):
        if  self.to_index!=self.index: self.to_index = self.index
    # move item/socket within 'technodes'/'inputs' up/down/to-index
    def upd_mv(self, context):
        if self.to_index==self.index: return  # Nothing to do
        inps=[i for i,skt in enumerate(self.inputs)if skt.bl_idname=="SktObjects"]
        idx, loc, len_oc =self.index, self.to_index, len(inps)
        if   loc < 0: loc = len_oc-1    # Wrap around from top to bottom
        elif loc > len_oc-1: loc = 0    # Wrap around from bottom to top
        if loc==idx: self.index = self.to_index = loc; return   # Nothing to do
        # gather all indices from the previous to the current/target zone-socket
        idxs = list(range(inps[idx-1]+1 if idx else 0, inps[idx]+1))# current
        idx_t= inps[loc-(loc<idx)]+(loc<idx) if loc else 0# target
        # we always start from the opposite end => rolls over, indices stay same
        for c in idxs: self.inputs.move(idxs[-(loc<idx)], idx_t)
        self.technodes.move(idx, loc)# move the collection item
        self.index = self.to_index  = loc # update index
    ### PROPERTIES
    id_string : bpy.props.StringProperty()# ui-list layout
    technodes : bpy.props.CollectionProperty(type=PG_Technode)
    index:  bpy.props.IntProperty(min=-1, default=-1, update=upd_idx)# currently
    # function-buttons (cheaply fake non-undo-able operators)
    to_index:bpy.props.IntProperty(min=-1, default=-1, update=upd_mv)# currently
    ptr_add: bpy.props.BoolProperty(update=upd_ptr_add)# links object in to list
    ptr_del: bpy.props.BoolProperty(update=upd_ptr_del)# remove object from list
    def draw_buttons(self, context, layout):
        row = layout.row(align=True)
        row.prop(self, "ptr_add", text="Add", icon="PLUS")
        row.prop(self, "ptr_del", text="Remove", icon="X")
        #layout.template_list("FILTERS_UL_simple", self.id_string,
        layout.template_list("TECHNODE_UL_simple", self.id_string,#self.path_from_id(),
            self, "technodes", self, "index", rows=2, maxrows=5)
        layout.prop(self, 'to_index', text="move to")# (=> update moves)
        ###########
        #return # drawn in the layout-list
        box = layout.column(align=True)
        # only draw for the active thing
        if  self.technodes and self.index!=-1:
            box.context_pointer_set('node', self)
            #op = box.operator("ifc2sim.update_technodes", text="", icon="FILE_REFRESH")
            #op.index, op.mode = self.index, "ALL"
            self.technodes[self.index].draw(self, box, self.index)
    ### we need: Add Socket? Tech-type, Zones?

# TODO: write context-node and simulation. Then: draw!
# => wir müssen eigentlich beim Geometrie-Export den VI-SUITE tree aufbauen.
#    Darin erfolgen die weiteren eingaben. DEFAULT= einfache Standard-Eingaben
#    UND: Was ist mit default-Materialien?
#### => Speicherkapazität Kühlung SimSIA
# Context-node als leeren Node, der nur das interface des echten Context-nodes streamt?
class NodeEnviContext(Node, IFC_SIM_Tree_Node):
    '''Sets up or creates the VI-Suite network and relays the interface elements'''
    bl_idname = "NodeEnviContext"
    bl_label = "Envi Export Context"
    vi_network : bpy.props.PointerProperty(type=bpy.types.NodeTree)
    ## PointerProperties cannot reference nodes. => fake it => unstable against renaming
    #vi_location : bpy.props.PointerProperty(type=bpy.types.Node)
    #envi_context : bpy.props.PointerProperty(type=bpy.types.Node)
    #envi_geometry : bpy.props.PointerProperty(type=bpy.types.Node)
    #envi_simulation : bpy.props.PointerProperty(type=bpy.types.Node)
    vi_location_name : bpy.props.StringProperty()
    @property
    def vi_location(self):
        return self.vi_network.nodes.get(self.vi_location_name, None)
    @vi_location.setter
    def vi_location(self, node): self.vi_location_name = node.name
    envi_context_name : bpy.props.StringProperty()
    @property
    def envi_context(self):
        return self.vi_network.nodes.get(self.envi_context_name, None)
    @envi_context.setter
    def envi_context(self, node): self.envi_context_name = node.name
    envi_geometry_name : bpy.props.StringProperty()
    @property
    def envi_geometry(self):
        return self.vi_network.nodes.get(self.envi_geometry_name, None)
    @envi_geometry.setter
    def envi_geometry(self, node): self.envi_geometry_name = node.name
    envi_simulation_name : bpy.props.StringProperty()
    @property
    def envi_simulation(self):
        return self.vi_network.nodes.get(self.envi_simulation_name, None)
    @envi_simulation.setter
    def envi_simulation(self, node): self.envi_simulation_name = node.name
    ##
    def init(self, context):
        self.fetch_references( context)# get the vi-network and its nodes
        #self.push_location => should we somehow enhance the location node?
        #self.push_geometry(self, context)# push the geometry-info into its node
    def fetch_references(self,_context=None):
        "Gathers the vi-network and its nodes. Where necessary, create it/them."
        for nt in bpy.data.node_groups:# fetch or set-up node tree
            if nt.bl_idname == 'ViN':
                self.vi_network = nt
                for nd in nt.nodes:
                    if nd.bl_idname=='No_En_Con': self.envi_context    = nd# context/export node
                    if nd.bl_idname=='No_Loc'   : self.vi_location     = nd
                    if nd.bl_idname=='No_En_Geo': self.envi_geometry   = nd
                    if nd.bl_idname=='No_En_Sim': self.envi_simulation = nd
                break
        else:
            self.vi_network = bpy.data.node_groups.new("VI-Network",'ViN')
        if not self.vi_location:
            self.vi_location    = self.vi_network.nodes.new('No_Loc'   )# 'VI Location'
            self.vi_location.loc= "1"
        if not self.envi_context:
            self.envi_context   = self.vi_network.nodes.new('No_En_Con')# 'EnVi Export'
        if not self.envi_geometry:
            self.envi_geometry  = self.vi_network.nodes.new('No_En_Geo')# 'EnVi Geometry'
        if not self.envi_simulation:
            self.envi_simulation= self.vi_network.nodes.new('No_En_Sim')# 'EnVi Simulation'
        ### setup links:
        skt_geo = self.envi_context.inputs[ 'Geometry in']
        skt_loc = self.envi_context.inputs[ 'Location in']
        skt_sim = self.envi_context.outputs['Context out']
        if not skt_geo.is_linked:
            self.vi_network.links.new(self.envi_geometry.outputs[0], skt_geo)
        if not skt_geo.is_linked:
            self.vi_network.links.new(self.vi_location.outputs[0], skt_loc)
        if not skt_sim.is_linked:
            self.vi_network.links.new(skt_sim, self.envi_simulation.inputs[0])
    # NOTE: The Geometry-Node has no properties nor id-properties necessary to relay.
    def push_geometry(self, context): pass# push the geometry-info into its node
    setup_context : bpy.props.BoolProperty(set=fetch_references,
        description="Fetch (or create) the VI-Network and its nodes.")
    def draw_buttons(self, context, layout):
        if not self.envi_context:
            b = layout.box()
            b.alert=True
            b.label(text="Context not set up!")
            b.label(text="Set it up with this button:")
            b.prop(self, 'setup_context', text="Fetch/Create context", toggle=True)
            return
        # location section
        layout.context_pointer_set("node", self.vi_location)
        if self.vi_location.loc!="1":
            b = layout.box()
            b.alert=True
            b.label(text="Please set the source to 'EPW'")
            self.vi_location.draw_buttons(context, b)
        else:
            layout.prop(self.vi_location, "weather", text="Weather File")
        # context section
        layout.context_pointer_set("node", self.envi_context)
        self.envi_context.draw_buttons(context, layout)
        # simulation section
        #if not self.envi_context.use_custom_color:
        layout.context_pointer_set("node", self.envi_simulation)
        self.envi_simulation.draw_buttons(context, layout)

# chart-node must feature:
#   0-2 input-sockets (results, objects-to-restrict-to)
#   multiple inputs: it's probably better supplied as collection, but about 12
#       fixed inputs would be enough. Unused inputs must be hidden.
#   scaling, axis: two Y-axes, maybe even multi-spined axes like here:
#       https://matplotlib.org/3.4.3/gallery/ticks_and_spines/multiple_yaxis_with_spines.html
#       and here: https://stackoverflow.com/questions/1982770/matplotlib-changing-the-color-of-an-axis
#   saving: ..
class NodeChart(Node, IFC_SIM_Tree_Node):
    '''Chart the input data'''
    bl_idname   = 'NodeChart'
    bl_label    = 'Chart'

    # collection of calls?
    #

# TODO: cache the readeso-handle?
esoreaders = {}# we store esoreaders in here
class NodeAnalyseZones(Node, IFC_SIM_Tree_Node):
    '''Reads results of simulated zones. From these, the zones are then sorted.'''
    bl_idname   = 'NodeAnalyseZones'
    bl_label    = 'Analyse Zones'
    bl_icon     = 'VIEWZOOM'
    # pointer to EnVi-Network NOTE:  if you pointer defines a specific subtype
    #   you cant change it from GUI, only from Python and only to said subtype.
    #envi_network : bpy.props.PointerProperty(type=bpy.types.NodeTree)
    id_string   : bpy.props.StringProperty()# for referencing purposes
    zones = None
    def copy(self, node):
        self.id_string = str(id(self))# setup id
    def free(self):
        if self.id_string in esoreaders: esoreaders.pop(self.id_string)
    def init(self, context):
        self.id_string = str(id(self))# setup id
        self.outputs.new('SktObjects',"Zones")
        #self.outputs.new('','sorted zones')
        # set up EnVi-Network pointer (If none exists, create one)
    @property
    def _esoreader_is_stale(self): # for internal use only ...
        return ( self.id_string not in esoreaders or
            esoreaders[self.id_string].fname!=self.eso_file )
    def get_esoreader(self):
        if self._esoreader_is_stale: # either no reader or to prevent stale data
            esoreaders[self.id_string] = readeso.ESOReader(self.eso_file)
            esoreaders[self.id_string].read()
        return esoreaders[self.id_string]
    #zones = list(set([x[0] for x in data.columns if x[0].startswith('EN_') and 'AIR' not in x[0]]))
    def upd_sort_by(self, context):
        def mean(arraylike):
            x=0
            for i,y in enumerate(arraylike):
                x = x + (y-x)/(i+1)
            return x
        if   self.sort_by=="MAXIMUM":    fnk=max
        elif self.sort_by=="MINIMUM":    fnk=min
        elif self.sort_by=="DEVIATION":
            def fnk(arraylike):
                x = mean(arraylike)
                return max(abs(y-x) for y in arraylike)
        elif self.sort_by=="MITTELWERT": fnk=mean
        #
        if   self.sort_by_inner=="T_AIR" :       fnk2=lambda z:z['temp_air']
        elif self.sort_by_inner=="Q_HEAT":       fnk2=lambda z:z['q_heat_sens']
        #elif self.sort_by_inner=="Q_COOL":       fnk2=lambda z:z['temp_air']
        elif self.sort_by_inner=="INFILTRATION": fnk2=lambda z:z['inf_ach']
        else: fnk2 = lambda z:z[self.sort_by_inner]
        #
        data = self.get_esoreader().dataframes['Hourly']
        zones = set()
        for x in data.columns:
            if x[0].startswith('EN_') and 'AIR' not in x[0]:
                zones.add(x[0])
        zobj = lambda z:context.blend_data.objects[z]
        zones = sorted(zones, key=lambda z:fnk(fnk2(data[z])))
        self.outputs[0].objects = map(zobj, zones)
    def itm_sort_by(self, context):
        return [(s.upper(),s,s)for s in "Maximum, Minimum, Deviation, Mittelwert".split(', ')]
    sort_by       : bpy.props.EnumProperty(items=itm_sort_by,       update=upd_sort_by)
    def itm_sort_by_inner(self, context):
        return [(s,_dct.get(s,s),_dct.get(s,s))for s in self._get_fields()[2]]
        #return [(s.upper(),s,s)for s in "T_Air, Q_heat, Infiltration".split(', ')]
    sort_by_inner : bpy.props.EnumProperty(items=itm_sort_by_inner, update=upd_sort_by)
    def _get_fields(self):# how do I cache the fields?
        reader = self.get_esoreader()
        #eso_tree = reader.tree
        #fields = { sensor for tn in eso_tree for kat in eso_tree[tn]
            #for sensor in eso_tree[tn][kat] }
        #dframes = {}# these must be the top-level keys and therefore a subset of
        ## set(['Daily', 'Each', 'Hourly', 'Monthly', 'RunPeriod', 'TimeStep'])
        #kategories = {}# mid-level keys
        #fields = {}# Note that not all fields may work with all kategories
        #for tn in eso_tree:
            #dframes.add(tn)
            #for kat in eso_tree[tn]:
                #kategories.add(kat)
                #for sensor in eso_tree[tn][kat]:
                    #fields.add(sensor)
        #return dframes, kategories, fields
        dframes = reader.dataframes.keys()
        dframe  = reader.dataframes[self.dfkey]
        zones   = {col[0] for col in dframe.columns
                if col[0].startswith('EN_') and 'AIR' not in col[0]}
        fields  = {fld for z in zones for fld in dframe[z]}
        return dframes, zones, fields
    dfkey = "Hourly"
    #dfkey = bpy.props.EnumProperty(items=lambda s,c:sorted(s._get_fields[0],
                                    #key =lambda k:k=='Hourly'))
    def get_file(self, value):  # value => bool, unbenötigt
        fpth = bpy.context.blend_data.filepath# this is the path to the blender file
        dpth = fpth[:-6]# this is the path to the simulationfolder
        if not os.path.isdir(dpth): self.errors[0]=True;return#Not yet simulated
        root, dirs, files = next(os.walk(dpth))#
        eso  = next((f for f in files if f.endswith('.eso')), None)
        if not eso:                 self.errors[1]=True;return# os.path.isfile(pth)
        self.eso_file = os.path.join(dpth, eso)
        #else: raise KeyError("Cant find simulation results. you may input the path yourself")
    def upd_eso_file(self, _context):
        eso_file = self.eso_file
        eso_file2= bpy.path.abspath(eso_file)
        if eso_file!=eso_file2:self.eso_file=eso_file2
        self.errors[2] = not os.path.isfile(eso_file)# We couldn't find the input path
    eso_file    : bpy.props.StringProperty(subtype="FILE_PATH", update=upd_eso_file)
    def upd_eso_file_fetch(self, context):
        if self.eso_file_fetch: self.eso_file_fetch=False;return
        self.get_file(context)
    eso_file_fetch : bpy.props.BoolProperty(set=get_file)#update=get_file)=>update braucht context, set nicht
    errors : bpy.props.BoolVectorProperty(size=4)# No folder, No file, no eso file, couldn't find the zones
    def draw_buttons(self, context, layout):
        row = layout.row(align=True)
        row.prop(self, 'eso_file_fetch', text="", icon='FORWARD', toggle=True)
        row.prop(self, 'eso_file', text="")
        #row.prop(self, 'eso_file_fetch', text="", icon='FILE_REFRESH', toggle=True)
        #row.operator("screen.userpref_show")
        #self.outputs[0]._draw(context, layout, self, "Objects")
        layout.prop(self, 'sort_by', text="order")
        layout.prop(self, 'sort_by_inner', text="quantity")
        if not self.eso_file:
            #layout.ui_units_x=10# sets a fixed layout size, so sub-elements may even extend beond nodes...
            #sublayout.ui_units_y=20# no effect in uppermost layout.
            b = layout.box()
            b.alert=True
            #b.label(text=("32x.16:1234567890"*11)[:int((self.width-32)*.16)])
            #b.label(text="{}".format(self.width))
            b.label(text="! No ESO-file (results) set !")
            b.label(text="After simulation the eso-file can be fetched with the")

        if self.eso_file and os.path.isfile(self.eso_file):
            if self.outputs[0].objects_collection:
                layout.template_list("POINTERS_UL_general", self.id_string,
                    self.outputs[0], "objects_collection",  self.outputs[0], "index",
                rows=2, maxrows=5)
        else:
            b = layout.box()
            b.alert=True
            b.label(text="! Eso-File not found !")
            b.label(text="Maybe you didn't simulate?")
            b.label(text="Maybe you renamed or moved the project?")
    #    row = layout.row(align=True)
    #    row.label("outputs")
    #    row.prop(self, 'ptr_add', text="", icon='PLUS', toggle=True)
    #    row.prop(self, 'ptr_del', text="", icon='X',    toggle=True)
    #



################################################################################

#==============================================================================#
# TODO: rewrite for a proper operator. See Log.txt for how messy engexport is.
# NOTE: bpy.ops.material.envi_node() should be run on the materials?
class NodeExportEnViGeometry(Node, IFC_SIM_Tree_Node):
    '''Node for exporting Zones to EnVi-Geometry'''
    bl_label = "Export EnVi Geometry"
    bl_icon  = "HOME"
    def init(self, context): # this allows to prepare the input objects first
        self.inputs.new( 'SktObjects', "Zones"  ).force_draw=-1
        self.inputs.new( 'SktObjects', "Shading").force_draw=-1

    # The geometry-export allows for this offset-vector (pregeo calls it in)
    geo_offset: bpy.props.FloatVectorProperty(name="Offset",
        description="Offset the exported geometry", subtype='TRANSLATION')
    def preexport( self, scene):pass# Only to satisfy the operators interface
    def postexport(self):       pass# as node.engexport expects these to exist
    def upd_setup_zones(self, context):
        if self.setup_zones: self.setup_zones=False; return
        zones = self.inputs["Zones"  ].objects
        shades= self.inputs["Shading"].objects
        warnings = []
        for obj in shades:
            obj.vi_params.vi_type   = '1' # we're an envi-object (shading)
            obj.vi_params.envi_type = '1' # we're shading
        for obj in zones:
            obj.vi_params.vi_type   = '1' # we're an envi-object (zone)
            obj.vi_params.envi_type = '0' # we're construction
        for obj in zones+shades:
            ### objects might be linked to multiple collections.
            outer_collections, inner_collections = [], []
            for col in obj.users_collection:
                [outer_collections, inner_collections][
                    col.vi_params.envi_zone == 1].append(col)
            # either all our collections are zones or we are unlinked. If assume
            # the former, we must already be set where we belong. (otherwise ...
            if not outer_collections: continue  # ... we got a problem anyways)
            # check the inner_collections for our name.
            for col in inner_collections:
                if col.name==obj.name: break    # found my wrapper-collection.
            else:# No envi-zone-collection defined for us
                col = bpy.data.collections.new(obj.name)
                if col.name != obj.name: # only explain if we've got a new name
                    print("Added collection '%s' for object '%s'"%(col.name, obj.name))
                col.objects.link(obj)
            for c in outer_collections:
                c.children.link( col)
                c.objects.unlink(obj)
    setup_zones: bpy.props.BoolProperty(update=upd_setup_zones, name="Setup Zones")

    def draw_buttons(self, context, layout):
        #=> this could/should be an operator with an undo function! Complexity!
        layout.prop(self, "setup_zones", toggle=True)
        # geo_offset
        layout.prop(self, 'geo_offset')
        # TODO: Maybe we better context_pointer_set to the geometry-export node
        #   instead of presenting this node as geometry-export node? Compare to
        #   the way the context-node does it by relaying the layout. This would
        #   mean that fetching the nodes/vi-network must be done earlier/here!
        layout.operator('node.engexport', text="Export Geometry")

class NodeLinkExportedZones(Node, IFC_SIM_Tree_Node):
    """Iterate over all EnVi-Networks and link up the zones with each other"""
    bl_label = "Link exported zones"
    bl_icon  = "LINKED"
    def upd_link_exported_zones(self, context):
        if self.link_exported_zones: self.link_exported_zones=False; return
        link_exported_zones()# see => !!! FUNCTIONS !!!
    link_exported_zones: bpy.props.BoolProperty(update = upd_link_exported_zones)
    def draw_buttons(self, context, layout):
        #=> this should be an operator with an undo function! Complexity!
        layout.prop(self, 'link_exported_zones', text="Link Exported Zones", toggle=True)

# TODO: We should combine this with the node for exporting zones. That would
#       allow faster calculation and we could just chain it to the export...
class NodeSetVolume(Node, IFC_SIM_Tree_Node):
    # currently we recalculate&set the volumes for all exported zones
    bl_label = "Reset Volumes"
    bl_icon  = "HOME"
    def upd_reset_volumes(self, context):
        if self.reset_volumes: self.reset_volumes=False;return
        set_correct_zone_volume()
    reset_volumes   : bpy.props.BoolProperty(update=upd_reset_volumes)

    def draw_buttons(self, context, layout):
        #=> this could/should be an operator with an undo function! Complexity!
        layout.prop(self, 'reset_volumes', text="Reset Volumes", toggle=True)

#==============================================================================#
# TODO: Rewrite for a proper operator. See Log.txt for how messy engexport is.
# NOTE: bpy.ops.material.envi_node() should be run on the materials?
# TODO: Combined node for preparing, exporting and linking up zones, and setting
#       their volumes. Optimize to allow faster calculation and chaining.
class NodeEnViGeometryExport(Node, IFC_SIM_Tree_Node):
    '''Node for exporting Zones to EnVi-Geometry'''
    bl_label = "Export EnVi Geometry"
    bl_icon  = "HOME"
    def init(self, context): # this allows to prepare the input objects first
        self.inputs.new( 'SktObjects', "Zones"  ).force_draw=-1
        self.inputs.new( 'SktObjects', "Shading").force_draw=-1

    # The geometry-export allows for this offset-vector (pregeo calls it in)
    geo_offset: bpy.props.FloatVectorProperty(name="Offset",
        description="Offset the exported geometry", subtype='TRANSLATION')
    def preexport( self, scene):pass# Only to satisfy the operators interface
    def postexport(self):       pass# as node.engexport expects these to exist
    def upd_setup_zones(self, context):
        if self.setup_zones: self.setup_zones=False; return
        zones = self.inputs["Zones"  ].objects
        shades= self.inputs["Shading"].objects
        warnings = []
        for obj in shades:
            obj.vi_params.vi_type   = '1' # we're an envi-object (shading)
            obj.vi_params.envi_type = '1' # we're shading
        for obj in zones:
            obj.vi_params.vi_type   = '1' # we're an envi-object (zone)
            obj.vi_params.envi_type = '0' # we're construction
        for obj in zones+shades:
            ### objects might be linked to multiple collections.
            outer_collections, inner_collections = [], []
            for col in obj.users_collection:
                [outer_collections, inner_collections][
                    col.vi_params.envi_zone == 1].append(col)
            # either all our collections are zones or we are unlinked. If assume
            # the former, we must already be set where we belong. (otherwise ...
            if not outer_collections: continue  # ... we got a problem anyways)
            # check the inner_collections for our name.
            for col in inner_collections:
                if col.name==obj.name: break    # found my wrapper-collection.
            else:# No envi-zone-collection defined for us
                col = bpy.data.collections.new(obj.name)
                if col.name != obj.name: # only explain if we've got a new name
                    print("Added collection '%s' for object '%s'"%(col.name, obj.name))
                col.objects.link(obj)
            for c in outer_collections:
                c.children.link( col)
                c.objects.unlink(obj)
    setup_zones: bpy.props.BoolProperty(update=upd_setup_zones, name="Setup Zones")
    #   Iterate over all EnVi-Networks and link up the zones with each other
    link_zones   : bpy.props.BoolProperty(set=lambda s,v:link_exported_zones()) # see => !!! FUNCTIONS !!!
    reset_volumes: bpy.props.BoolProperty(set=lambda s,v:set_correct_zone_volume())

    def draw_buttons(self, context, layout):
        #=> this could/should be an operator with an undo function! Complexity!
        layout.prop(self, "setup_zones", toggle=True)
        # TODO: Maybe we better context_pointer_set to the geometry-export node
        #   instead of presenting this node as geometry-export node? Compare to
        #   the way the context-node does it by relaying the layout. This would
        #   mean that fetching the nodes/vi-network must be done earlier/here!
        layout.prop(self, 'geo_offset')
        layout.operator(  'node.engexport',   text="Export Geometry")
        #=> this could/should be an operator with an undo function! Complexity!
        layout.prop(self, 'link_zones', text="Link Exported Zones", toggle=True)
        layout.prop(self, 'reset_volumes',    text="Reset Volumes", toggle=True)
#==============================================================================#

# NOTE: we have no direct way of getting the imported objects (yet)
#       Therefore we define an operator that keeps track. Clearing optional.
class NodeImportIFC(Node, IFC_SIM_Tree_Node):
    """Iterate over all EnVi-Networks and link up the zones with each other"""
    bl_label = "IFC-Import"
    bl_icon  = "IMPORT"
    def init(self, context):
        self.outputs.new('SktObjects', "All Objects")
    def draw_buttons(self, context, layout):
        #row = layout.row()# we override context.selected_objects
        #row.context_pointer_set("selected_objects", list(context.view_layer.objects))
        #row.operator("object.delete", text="Delete Objects")
        #bpy.ops.object.delete({"selected_objects": objs})
        #=> this should be an operator with an undo function! Complexity!
        layout.operator('ifc2sim.import_ifc')

################################################################################
### Node to define a construction TODO: materials, constructions
class NodeConstruction(Node, IFC_SIM_Tree_Node):
    #bl_idname= "ifc2sim_construction"#
    bl_label = "Construction" # currently we only allow simple assumptions
    bl_icon  = "MATERIAL" # use an icon for this node?
    def init(self, context):
        self.outputs.new('NodeSocketVirtual', "Construction")# one output
        self.inputs.new( 'NodeSocketVirtual', "Material/Construction")# any number of materials or construction
        "outside material"
    material : bpy.props.PointerProperty(type=bpy.types.Material)
    def upd_set_material(self, context):
        if  self.set_materials:
            self.set_materials=False; return
        mat = self.material
        if not mat: return
        vip = self.material.vi_params#
        if not vip.envi_nodes: # we must create the node tree
            pass
            vip.envi_nodes= bpy.data.node_groups.new(mat.name, 'EnViMatN')
            vip.envi_nodes.nodes.new('No_En_Mat_Con').active = True
            vip.envi_nodes['envi_con_type']= 'None'
            vip.envi_nodes['enmatparams']  = {'airflow':0, 'boundary':0, 'tm':0}
        elif mat.name != vip.envi_nodes.name \
            and vip.envi_nodes.name in bpy.data.materials:# make a copy
            vip.envi_nodes = vip.envi_nodes.copy()
            vip.envi_nodes.name = mat.name
        cons= filter(lambda n:n.bl_idname=='No_En_Mat_Con',vip.envi_nodes.nodes)#
        cons= sorted(cons, key=lambda n:n.active)
        if cons: con = cons[-1]
        else:
            con = vip.envi_nodes.nodes.new('No_En_Mat_Con')
            con.active = True
        # Hier sollten jetzt die materialien aufgebaut werden...
    set_material : bpy.props.BoolProperty(update=upd_set_material)

    def draw_buttons(self, context, layout):
        #=> this should be an operator with an undo function! Complexity!
        #op = layout.operator("IFC2SIM.set_materials")
        layout.prop(self, 'set_materials', text="Set Materials", toggle=True)
        for c in ("Ceiling","Floor","North","East","South","West"):
            layout.prop(self, 'material_'+c.lower(), text=c)

# wie halte ich referenz? sont ewiges filtern...
###class NodeBSDummy(Node, IFC_SIM_Tree_Node):
    ###bl_label = "NodeBSDummy"
    ###def get_node(self): pass
        ####ln  = self.outputs[0].links[0]
        ####nd  = ln.to_node
        ####skt = ln.to_socket
    ###def draw_buttons(self, context, layout):
        ###for i,inp in enumerate(self.node.inputs):

            ###if inp.is_linked:
                ###node = inp.links[0].from_node
                ###node.draw_buttons(self, context, layout)

# END ##########################################################################
###                            !!! OPERATORS !!!                             ###
#######################################80############################### BEGIN #
### first: the function TODO: rewrite to multiple actives and to lowlevel
#   FIXME: Should we A) join the cutters togethter first and join that in to obj
#       or should we B) test each opening if it had cut, removing it if unneeded
#       Each excludes the other. A) is far simpler. How would B check if cutters
#       are still needed on other zones?
## something doesn't yet work on the first run, but does on the second...
## seems to be an stability-error within intersect. (ignored for now).
def cut_blocks_new(zones, objs, collection, material=None, vg_name=None, keep=(True,False), remove_vg=True):
    """ Cut in faces, optionally add them to a material, a vertex group or both """
    # NOTE: we'll assume, we're in objectmode
    # back up previous selection & active object # breaks if object gets deleted
    _selected   = bpy.context.selected_objects   # (e.g. if they are joined into
    _active     = bpy.context.view_layer.objects.active # any other mesh )
    # deselect all #operator: bpy.ops.object.select_all(action='DESELECT')
    for obj in _selected: obj.select_set(False)
    # copy objects, if requested:
    if keep[0]:
        objs = [obj.copy()for obj in objs]  # get a shallow/linked copy first
        for obj in objs:
            obj.data = obj.data.copy()      # make it stand-alone (copy mesh)
            collection.objects.link(obj)    # link it, so it won't get lost
    if keep[1]:
        zones = [obj.copy()for obj in zones]# get a shallow/linked copy first
        for obj in zones:
            obj.data = obj.data.copy()      # make it stand-alone (copy mesh)
            collection.objects.link(obj)    # link it, so it won't get lost
    # bundle up the cutters into one object:
    bpy.context.view_layer.objects.active = objs[0]
    if len(objs)>1: # only one cutter, no need to join.
        for obj in objs: obj.select_set(True)
        bpy.ops.object.join()
    cutter = bpy.context.view_layer.objects.active
    cutter.select_set(False)# deselect the cutter-object so it stays back
    # set mesh-selection mode to 'edge' # NOTE lowlevel, works in objectmode too
    bpy.context.tool_settings.mesh_select_mode[:] = [False,True,False]
    # loop over zones:
    for zone in zones:
        print("acting on zone ",zone.name)
        ### prepare room
        # select & activate; add vertex_group; add vertices to vertex_group;;
        zone.select_set(True); bpy.context.view_layer.objects.active = zone
        VG = zone.vertex_groups.new(name="cut_blocks__macro__temporary_vertex_group")
        indices = range(len(zone.data.vertices))# we want the indices for ALL vertices
        VG.add(indices, 1, "ADD")# adds the vertices to vertexgroup#=>objectmode
        # apparantly one may stack bmesh.from_mesh calls to join? => doesn't respect transforms
        # FIXME: if we use only one cutter, the VG could belong to it instead.

        ### join the cutter-object into this object, so we may use intersect.
        # TODO-Dummy first join, then select all vertices of pgroup
        # copy&link&select cutter. NOTE: since we join the cutter into the zone,
        #   a thin copy is enough. The other way round, we'd need a full copy.
        c_ob=cutter.copy(); collection.objects.link(c_ob); c_ob.select_set(True)
        print(len(zone.data.vertices), len(c_ob.data.vertices))
        # join the cutter into the zone
        bpy.ops.object.join()
        # select only vertices that are part of the vertex_group
        for v in zone.data.vertices:
            v.select = any(VG.index==g.group for g in v.groups)
        print(sum(v.select for v in zone.data.vertices), sum(not v.select for v in zone.data.vertices))

        ### now for the main task: cut the elements in
        # enter editmode; cut; enter object mode to update vertex-info
        bpy.ops.object.mode_set(mode='EDIT')
        ## set selection style to edge (needs editmode)# NOTE Done before loop
        #bpy.ops.mesh.select_mode(use_extend=False, use_expand=False, type='EDGE')
        # intersect. threshold=0 to prevent welding the pieces ## other ways?
        bpy.ops.mesh.intersect(mode='SELECT_UNSELECT', separate_mode='CUT', threshold=0.00)
        # check if we actually cut anything
        bpy.ops.object.mode_set(mode='OBJECT') # vertex-info actualisieren
        bpy.ops.object.mode_set(mode='EDIT') ## geht vermutlich mit dem dependency-graph
        #https://blender.stackexchange.com/questions/146559/how-do-i-get-a-mesh-data-block-with-modifiers-and-shape-keys-applied-in-blender
        #
        # apparantly we test for each selected vertex wether the VertexGroup is in it?
        _cut = any(VG.index not in [g.group for g in v.groups]
            for v in zone.data.vertices if v.select)
        ### if we did cut, we want to add material or vertexgroup, etc.
        if _cut: # abort if we did not cut
            # select inside faces (there's probably no lowlevel call for this)
            bpy.ops.mesh.loop_to_region(select_bigger=False)
            # add material; add vertex_group # FIXME: how to lowlevel ?
            if material==-1:
                zone.active_material_index = len(zone.material_slots)-1
                bpy.ops.object.material_slot_assign()
            elif material:
                bpy.ops.object.material_slot_add()
                bpy.ops.object.material_slot_assign()
                matSlot = zone.material_slots[zone.active_material_index]
                if type(material)==bpy.types.Material:
                    mat, material = material, material.name
                else:# name
                    mat = bpy.data.materials.get(material)
                if mat==None:
                    mat = bpy.data.materials.new(material)
                matSlot.material = mat
            if vg_name:
                ## add the selected vertices to new vertexgroup # => objectmode!
                #zone.vertex_groups.new(name=vg_name).add( [v.index
                #    for v in zone.data.vertices if v.select), 1, "ADD")
                ## this works in editmode, but needs operators
                bpy.ops.object.vertex_group_add()
                zone.vertex_groups[zone.vertex_groups.active_index].name=vg_name
                bpy.ops.object.vertex_group_assign()

        ### clean-up job.
        bpy.ops.mesh.select_all(action='DESELECT')  ## remove all that came from
        zone.vertex_groups.active_index = VG.index  #  the cutter-objects
        bpy.ops.object.vertex_group_select()        # \
        bpy.ops.mesh.select_linked(delimit=set())   #  } select not-in-vertex_group
        bpy.ops.mesh.select_all(action='INVERT')    # /
        bpy.ops.mesh.delete(type='VERT')            #  deletion needs editmode !!!
        if remove_vg: zone.vertex_groups.remove(VG) ## remove vertex-group
        bpy.ops.object.mode_set(mode="OBJECT")      ## restore object mode
        zone.select_set(False) ## prevent erraneous selection in join processes

    ### clean-up / restore previous state ## skipped. (many try-excepts)
    ## try restoring active object. if unavailable, unset active object.
    #try:    bpy.context.view_layer.objects.active = active
    #except ReferenceError as RE: bpy.context.view_layer.objects.active = None
    ## try restoring selection
    return zones, cutter # => we return the cutter, just in case


def cut_blocks(active=None, objs=None, material=None, vertex_group=None, keep=(True,False), remove_vg=True):
    """ Cut in Faces, optionally add them to a material, a vertex group or both
    :active: object to cut into (e.g. room), defaults to active object
    :objs: list of cutters (e.g. windows), defaults to selection
    :material: add new faces to the material named here. If it doesn't exist, create it.
    :vertex_group: like :material: but for vertex groups
    :keep: 2*bool - keep copies of the provided cutters and active, defaults to (True,False)
    :remove_vg: remove the base-vertexgroup marking the e.g. room, defaults to true
    """
    # prepare context
    ctxt = bpy.context.copy()
    ctxt['active_object'] = ctxt['scene'].objects[0]
    try:
        bpy.ops.object.mode_set(ctxt, mode="OBJECT")
    except:
        pass # allready in objectmode...

    # get objects
    if not active: active = bpy.context.view_layer.objects.active
    if not objs: objs = [obj for obj in bpy.context.scene.objects if obj.select_get() and obj!=active]

    bpy.ops.object.select_all(action='DESELECT')
    if keep[0]:
        for obj in objs: obj.select_set(True)
        bpy.ops.object.duplicate(linked=False, mode='TRANSLATION')
        objs = [obj for obj in bpy.context.scene.objects if obj.select_get()]
    bpy.ops.object.select_all(action='DESELECT')
    if keep[1]:
        active.select_set(True)
        bpy.ops.object.duplicate(linked=False, mode='TRANSLATION')
        active = [obj for obj in bpy.context.scene.objects if obj.select_get()][0]
    bpy.ops.object.select_all(action='DESELECT')

    _active = bpy.context.view_layer.objects.active
    bpy.context.view_layer.objects.active = active

    # prepare room
    active.select_set(True)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.object.vertex_group_add()
    vg_base=active.vertex_groups[active.vertex_groups.active_index]
    vg_base.name="cut_blocks__macro__temporary_vertex_group"
    bpy.ops.object.vertex_group_assign()
    bpy.ops.object.mode_set(mode="OBJECT")

    for obj in objs: obj.select_set(True)
    bpy.ops.object.join()
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.vertex_group_select()


    bpy.ops.mesh.select_mode(use_extend=False, use_expand=False, type='EDGE')
    #bpy.ops.mesh.intersect(mode='SELECT_UNSELECT', separate_mode='CUT', threshold=1e-06)
    # 0-threshold. otherwise it could weld together both pieces ## Anderer weg?
    bpy.ops.mesh.intersect(mode='SELECT_UNSELECT', separate_mode='CUT', threshold=0.00)

    # check if we actually did cut anything?
    bpy.ops.object.mode_set(mode='OBJECT') # vertex-info actualisieren
    bpy.ops.object.mode_set(mode='EDIT')
    _cut = any(vg_base.index not in [g.group for g in v.groups]
        for v in active.data.vertices if v.select)

    if _cut: # abort if we did not cut
        bpy.ops.mesh.loop_to_region(select_bigger=False)

        if material==-1:
            active.active_material_index = len(active.material_slots)-1
            bpy.ops.object.material_slot_assign()
        elif material:
            bpy.ops.object.material_slot_add()
            bpy.ops.object.material_slot_assign()
            matSlot = active.material_slots[active.active_material_index]
            mat = bpy.data.materials.get(material)
            if mat==None:
                mat = bpy.data.materials.new(material)
            matSlot.material = mat
        if vertex_group:
            bpy.ops.object.vertex_group_add()
            group=active.vertex_groups[active.vertex_groups.active_index]
            group.name=vertex_group
            bpy.ops.object.vertex_group_assign()

    # cleanup
    bpy.ops.mesh.select_all(action='DESELECT')
    active.vertex_groups.active_index = vg_base.index
    bpy.ops.object.vertex_group_select()
    bpy.ops.mesh.select_linked(delimit=set())
    bpy.ops.mesh.select_all(action='INVERT')
    bpy.ops.mesh.delete(type='VERT')

    # ungenutzte material-slots entfernen
    #active.material_slots.remove ...

    if remove_vg:
        #bpy.ops.object.vertex_group_remove(all=False, all_unlocked=False)
        active.vertex_groups.remove(vg_base)

    bpy.ops.object.mode_set(mode="OBJECT")
    return active


# NOTE: Operators can NOT use PointerProperties or PointerGroups to Objects, etc
#   The ValueError reads 'bpy_struct "..." doesn't support datablock properties'
# NOTE: Also layout.operator returns just a Properties-bundle, skipping PyProps
# FIXME: Therefore we could store the names in a context_pointer_set, but since
#   the context already contains the spawning node as 'node', we may use that.
# TODO: How to make a non-node operator do the same?
class IFC2SIM_OT_cut_openings(bpy.types.Operator):
    """Operator to cut window-/door-openings into zones"""
    bl_idname = "ifc2sim.cut_openings"
    bl_label  = "Cut windows/doors into zones"
    bl_options= {"UNDO","INTERNAL","REGISTER"}

    @classmethod # FIXME: it'd be sensible to make this available in 3D somehow.
    def poll(cls, context):
        return (context.space_data.type == 'NODE_EDITOR'
            and context.space_data.tree_type == IFC_SIM_Tree.bl_idname)

    ### PROPERTIES # not yet used -- remove them?
    # we should use copies, as the originals use different names, etc
    use_copies: bpy.props.BoolProperty(default=True)
    # if we find objects at our targets, they're old => delete them
    remove_old: bpy.props.BoolProperty(default=True)

    def execute(self, context):
        node    = context.node
        doors   = node.inputs["Doors"   ].objects
        windows = node.inputs["Windows" ].objects
        zones   = node.inputs["Zones"   ].objects
        mat_door    = node.material_door
        mat_window  = node.material_window
        ### capsulated away to ensure access from outside
        state, zones= self.execute_internal(zones, doors, windows,
                                            mat_door, mat_window)
        if state=="CANCELLED": return {"CANCELLED"}
        ### return the zones after modifications are done.
        skt = node.outputs.get("Zones")     # get output socket, if available
        if not skt: skt = node.outputs.new('SktObjects', "Zones")# add if needed
        skt.objects = zones                 # push the modified zones.
        return {state}
    # this main part is put into a static method to allow *some* access from
    @staticmethod   # outside. Really to skip a op-layer for ifc2sim.setup_tree
    def execute_internal(zones, doors, windows, mat_door, mat_window):
        if not zones: return 'CANCELLED', []
        collection = zones[0].users_collection[0]# TODO: collection-input
        ### cut doors first and windows just after. That way if someone defined
        #   the glass of their door as window, we'll add it as such by default.
        if doors:
            print("CUTTING DOORS #############################################")
            zones, cutter = cut_blocks_new(zones, doors,   collection,
                                           mat_door)
            bpy.data.objects.remove(cutter) # cleanup-job: delete the cutters
        if windows:
            print("CUTTING WINDOWS ###########################################")
            zones, cutter = cut_blocks_new(zones, windows, collection,
                                           mat_window)
            bpy.data.objects.remove(cutter) # cleanup-job: delete the cutters
        return 'FINISHED', zones


# NOTE: I stole the contents of the bim.load_project operator, because calling
#       it by INVOKE_DEFAULT returns the invokes result, so we continue to early
#       as that would be RUNNING_MODAL. Instead we turn ourselves to a bloated
#       copy, as to provide the gui
# TODO: Maybe split the objects into their zone-types directly?
class IFC2SIM_OT_import_ifc(bpy.types.Operator):
    """Imports IFC-File and, if provided an output, returns new objects"""
    bl_idname = "ifc2sim.import_ifc"
    bl_label  = "Import IFC"
    bl_options= {"UNDO","REGISTER"}
    bl_description = "Load/Import IFC project"
    filepath: bpy.props.StringProperty(subtype="FILE_PATH")
    filter_glob: bpy.props.StringProperty(default="*.ifc;*.ifczip;*.ifcxml", options={"HIDDEN"})
    is_advanced: bpy.props.BoolProperty(name="Enable Advanced Mode", default=False)

    @classmethod # FIXME: it'd be sensible to make this available in 3D somehow.
    def poll(cls, context): # The main part doesn't depend on there being a node
        return (context.space_data.type == 'NODE_EDITOR'
            and context.space_data.tree_type == IFC_SIM_Tree.bl_idname)

    def execute(self, context):
        VL = context.view_layer
        objects = list(VL.objects)
        #bpy.ops.bim.load_project("INVOKE_DEFAULT")
        fn =bpy.ops.bim.load_project('EXEC_DEFAULT', filepath=self.filepath,
            filter_glob=self.filter_glob, is_advanced=self.is_advanced)
        if hasattr(self, "node"):       node = self.node
        elif hasattr(context, "node"):  node = context.node
        else: return {"FINISHED"}
        # we got a node from somewhere? return the objects
        all_objs = node.outputs.find("All Objects")
        if all_objs==-1:
            all_objs = node.outputs.new('SktObjects', "All Objects")
        else: all_objs=node.outputs[all_objs]
        all_objs.objects=[o for o in VL.objects if o not in objects]
        return {"FINISHED"}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        self.node = context.node if hasattr(context, "node") else None
        return {"RUNNING_MODAL"}

# TODO: How to make a non-node operator do the same?
class IFC2SIM_OT_collection_items(bpy.types.Operator):
    """This operator adds/removes/sets-to the provided objects/collections"""
    bl_idname = "ifc2sim.collection_items"
    bl_label  = "Modify objects/children of collection"
    bl_options= {"UNDO","INTERNAL","REGISTER"}

    @classmethod # FIXME: it'd be sensible to make this available in 3D somehow.
    def poll(cls, context):
        return (context.space_data.type == 'NODE_EDITOR'
            and context.space_data.tree_type == IFC_SIM_Tree.bl_idname)

    ### PROPERTIES
    single:bpy.props.BoolProperty(description="If true, update only current")
    mode : bpy.props.EnumProperty(items=[
        ("ADD",     "Add",      "Add input to collections children/objects" ),
        ("REMOVE",  "Remove",   "Remove from collection children/objects"   ),
        ("SET",     "Set",      "Set collection children/objects to input"  ),
        #("ADD_ALL",     "Add",      "For each input Add input to collections children/objects" ),
        #("REMOVE_ALL",  "Remove",   "For each input Remove from collection children/objects"   ),
        #("SET_ALL",     "Set",      "For each input Set collection children/objects to input"  )
        ])

    def execute(self, context):
        node = context.node
        mode = self.mode#.split("_")
        #single=len(mode)-1
        #mode = mode[0]
        for i,(skt_a,skt_b)in enumerate(zip(node.inputs[:-1], node.inputs[1:])):
            if node.single and not i==(node.index//2)*2: continue
            #if self.single and not i==(node.index//2)*2: continue
            #if single and not i==(node.index//2)*2: continue
            elif i%2: continue
            if skt_a.identifier.startswith("... to those Collections"):
                skt_a, skt_b = skt_b, skt_a
            collections = skt_b.collections
            links = skt_a.links
            for lnk in links:
                skt = lnk.from_socket
                if skt.bl_idname=="SktCollections":
                    colls= skt.collections
                    if   mode=='ADD':
                        for col in collections:
                            ch = col.children
                            for c in colls:
                                if not c.name in ch: ch.link(c)
                    elif mode=='SET':
                        for col in collections:
                            ch = col.children
                            for c in set(ch): ch.unlink(c)
                            for c in colls  : ch.link(  c)
                    elif mode=='REMOVE':
                        for col in collections:
                            ch = col.children
                            for c in colls:
                                if c.name in ch: ch.unlink(c)
                if lnk.from_socket.bl_idname=="SktObjects":
                    objs = skt.objects
                    if   mode=='ADD':
                        for col in collections:
                            ch = col.objects
                            for c in objs:
                                if not c.name in ch: ch.link(c)
                    elif mode=='SET':
                        for col in collections:
                            ch = col.objects
                            for c in set(ch): ch.unlink(c)
                            for c in objs   : ch.link(  c)
                    elif mode=='REMOVE':
                        for col in collections:
                            ch = col.objects
                            for c in objs:
                                if c.name in ch: ch.unlink(c)
        return {'FINISHED'}

# TODO: rebuilding the 'update'-thingies. execute 'handle'-thingy
class IFC2SIM_OT_update_technodes(bpy.types.Operator):
    """This operator adds/removes Tech-nodes and updates the EnVi-Network from them"""
    bl_idname = "ifc2sim.update_technodes"
    bl_label  = "Update"
    bl_options= {"UNDO","INTERNAL","REGISTER"}

    @classmethod
    def poll(cls, context):
        return (context.space_data.type == 'NODE_EDITOR'
            and context.space_data.tree_type == IFC_SIM_Tree.bl_idname)

    ### PROPERTIES
    index: bpy.props.IntProperty(description="Index of the SktObjects to work on")
    mode : bpy.props.EnumProperty(items=[(s,s,"Add/Remove/Update from "+s+"-node")
        for s in("HVAC", "Infiltration", "Equipment", "Occupancy" ) ]+[
        ("ALL","ALL","Add/Remove/Update from all tech-nodes")])

    def execute(self, context):
        node, mode, index = context.node, self.mode, self.index
        self.node = node
        sockets= [ins for ins in enumerate(node.inputs)if ins[1].bl_idname=="SktObjects"]
        if index>=len(sockets):
            self.report({'WARNING'},
                "Missmatch between number of SktObjects-sockets and items")
            return {'CANCELLED'}
        # const: index of previous SktObjects-socket to search for items between
        idx_a = sockets[index-1][0] if index>0 else 0
        skt_z = sockets[index][1]# der zonen-socket der betrachtet wird
        identifier = sockets[index][1].identifier
        for i, key, bl_idname in (
            (0,'HVAC','So_En_Net_Hvac'), (1,'Infiltration','So_En_Net_Inf'),
            (2,'Equipment','So_En_Net_Eq'), (3,'Occupancy','So_En_Net_Occ'), ):
            if self.mode in {'ALL', key}:
                item  = node.technodes[index];in_use=item.technodes[i]# in_use?
                ### while the docstring claims to use the identifier it uses the name
                #idx_b = node.inputs.find(identifier)# recalc up-to-date index
                ### thats why we loop directly
                for ioff, skt in enumerate(node.inputs[idx_a:]):
                    if skt.bl_idname==bl_idname:exists=True;break# we got one
                    if skt==skt_z:exists=False;break# we're sure there is none
                else: raise IndexError('The socket got lost on the way?')       # should be impossible
                idx_b = idx_a + ioff + 1 # recalc up-to-date index
                if  (exists,in_use)==(True,False): self.kill_input(skt)# remove
                elif(exists,in_use)==(False,True): # add new input
                    skt  = node.inputs.new(bl_idname, key)
                    idx_c= len(node.inputs)-1; node.inputs.move(idx_c, idx_b-1)
                    # add & setup node based on educated guess and then link it
                    self.add_input(skt)# to the socket (gets node-type from id)
                ## update from Socket/Node (unconnected Socket doesn't update)
                ## pass
                # finally set the mode, as to know we did the update already and
                item.technodes[4+i]=item.technodes[i]# don't need highlighting
        return {'FINISHED'}
    # need self for applying educated guess.
    def add_input(self, skt):
        pass
        d=dict([('So_En_Net_Hvac',('No_En_Net_Hvac',self.setup_HVAC)),
                ('So_En_Net_Inf', ('No_En_Net_Inf', self.setup_Infiltration)),
                ('So_En_Net_Eq',  ('No_En_Net_Eq',  self.setup_Equipment)),
                ('So_En_Net_Occ', ('No_En_Net_Occ', self.setup_Occupancy))      ])
        if skt.bl_idname in d:
            ## add node and link
            node = skt.node.id_data.nodes.new(d[skt.bl_idname][0])
            node.id_data.links.new(node.outputs[0], skt)
            node.hide=True;node.location=self.node.location;node.location.x-=200
            ## add educated guess
            #d[ skt.bl_idname ]( node, zones)
    @staticmethod # just a utility, no need for self.
    def kill_input(skt):
        for lnk in skt.links:
            skt_b = lnk.from_socket
            if len(skt_b.links)==1: skt.node.id_data.nodes.remove(skt_b.node)   # kill node if not needed
        skt.node.inputs.remove(skt)                                             # kills the socket
    # TODO: actually implement this ...
    def setup_HVAC(         self, node, zones): pass
    def setup_Infiltration( self, node, zones): pass
    def setup_Equipment(    self, node, zones): pass
    def setup_Occupancy(    self, node, zones): pass
    # TODO: actually implement this ...
    # use a staticmethod because we dont use the operator
    @staticmethod
    def update_HVAC(node, mode, index, idx_a, idx_b):
        for skt in node.inputs[idx_a:idx_b]:
            if skt.bl_idname=='So_En_Net_Hvac': break
        else:
            skt=node.inputs.new('So_En_Net_Hvac', "HVAC")
            idx_c=len(node.inputs)-1; node.inputs.move(idx_c, idx_b)
            # add & setup node based on educated guess and link it to the socket
            # pass
        # update from Socket/Node (unconnected Socket doesn't update)
        # pass
        # finally set the mode, as to know we did the update already and need
        item=node.technodes[index];item.technodes[4]=item.technodes[0]# no highlighting
    @staticmethod
    def update_Infiltration(node, mode, index, idx_a, idx_b):
        for skt in node.inputs[idx_a:idx_b]:
            if skt.bl_idname=='So_En_Net_Inf': break
        else:
            skt=node.inputs.new('So_En_Net_Inf', "Infiltration")
            idx_c=len(node.inputs)-1; node.inputs.move(idx_c, idx_b)
            # add & setup node based on educated guess and link it to the socket
            # pass
        # update from Socket/Node (unconnected Socket doesn't update)
        # pass
        # finally set the mode, as to know we did the update already and need
        item=node.technodes[index];item.technodes[5]=item.technodes[1]# no highlighting
    @staticmethod
    def update_Equipment(node, mode, index, idx_a, idx_b):
        for skt in node.inputs[idx_a:idx_b]:
            if skt.bl_idname=='So_En_Net_Eq': break
        else:
            skt=node.inputs.new('So_En_Net_Eq', 'Equipment')
            idx_c=len(node.inputs)-1; node.inputs.move(idx_c, idx_b)
            # add & setup node based on educated guess and link it to the socket
            # pass
        # update from Socket/Node (unconnected Socket doesn't update)
        # pass
        # finally set the mode, as to know we did the update already and need
        item=node.technodes[index];item.technodes[6]=item.technodes[2]# no highlighting
    @staticmethod
    def update_Occupancy(node, mode, index, idx_a, idx_b):
        for skt in node.inputs[idx_a:idx_b]:
            if skt.bl_idname=='So_En_Net_Occ': break
        else:
            skt=node.inputs.new('So_En_Net_Occ', "Occupancy")
            idx_c=len(node.inputs)-1; node.inputs.move(idx_c, idx_b)
            # add & setup node based on educated guess and link it to the socket
            # pass
        # update from Socket/Node (unconnected Socket doesn't update)
        # pass
        # finally set the mode, as to know we did the update already and need
        item=node.technodes[index];item.technodes[7]=item.technodes[3]# no highlighting
    @staticmethod
    def get_zones(node, idx_b):
        # our objects might either be unexported zones, or their exported copies
        # if exported, their only collection is of same name (EN_...) and direct
        #   child to the collection 'EnVi Geometry'
        # if not, we expect an exported collection for each of these collections
        #   here. NOTE: We skip non-exported collections in updating
        #   Should we report erraneous collections and-or objects?
        collections = set()
        zones = set(node.inputs[idx_b].objects)# first gather our objects
        for zone in zones: # then extend with their collections
            collections |= zone.users_collection
        missing_collections = set() # we store the missing zones for detail
        envi_zone_names = set()
        envi_exported = bpy.data.collections['EnVi Geometry'].children
        for coll in collections:
            c_name = 'EN_'+coll.name.upper().replace('-', '_').replace('/', '_')
            if coll.name in envi_exported: envi_zone_names.add(coll.name)
            elif c_name  in envi_exported: envi_zone_names.add(c_name)
            else: missing_collections.add(coll.name)
        ### we got the collections (or rather their names), now gather the nodes
        # it's unclear, wether we should support repeated zone nodes? probably
        # not, but sorting for the leftmost instance is stupid as well...
        zone_nodes = [n for n in node.envi_network.nodes
                        if  n.bl_idname=="No_En_Net_Zone"
                        and n.zone in envi_zone_names]
        missing_nodes = envi_zone_names - {n.zone for n in zone_nodes}#=> detail
        return zone_nodes, envi_zone_names, missing_nodes, missing_collections
    # now, the actual update needs to go to the node and see their inputs.
    # then in must check their type and update accordingly...
    @staticmethod
    def handle_tech_in_EnViN(node, mode, index, idx_a, idx_b):
        """Create/Update/Clear tech-nodes in the EnVi-Network"""
        ### the needed information of said mode:
        #   identifier      SOCKET:bl_idname    NODE:bl_idname   /settings
        __data = [('HVAC',  'So_En_Net_Hvac',   'No_En_Net_Hvac',
                {'envi_hvact',  'envi_hvacht',  'envi_hvacct',
                'envi_hvachlt', 'envi_hvachaf', 'envi_hvacshc',
                'envi_hvacclt', 'envi_hvaccaf', 'envi_hvacscc',
                'envi_hvacoam', 'envi_hvacfrp', 'envi_hvacfrzfa',
                'envi_hvacfrz', 'envi_hvacfach','envi_hvachr',
                'envi_hvachre', 'h', 'c',
                'acttype', 'envi_heat', 'envi_htsp',
                'envi_cool', 'envi_ctsp'}                           ),
            ('Infiltration','So_En_Net_Inf',    'No_En_Net_Inf',
                {'envi_inftype','envi_inflevel'}                    ),
            ('Occupancy',   'So_En_Net_Occ',    'No_En_Net_Occ',
                {'envi_occwatts','envi_weff',   'envi_airv',
                'envi_cloth',   'envi_occtype', 'envi_occsmax',
                'envi_comfort', 'envi_co2'}                         ),
            ('Equipment',   'So_En_Net_Eq',     'No_En_Net_Eq',
                {'envi_equiptype','envi_equipmax'}                  )]
        envi_network = node.envi_network
        zones, _names, _zones_amiss, _colls_amiss = self.get_zones(node, idx_b)
        if not zones or not envi_network: return # exit-early
        zones = sorted(zones, key=lambda z:z.location)
        X0,Y0 = zones[0].location
        # setup HVAC
        for idx_m, (MODE, skt_bl_idname, tnd_bl_idname, settings) in enumerate(__data):
            if mode not in {'ALL', MODE}: continue # early-exit
            add = node.technodes[index].technodes[idx_m]# True/False=>add/remove
            ### get tech-nodes if existent, as well as zone-nodes for speed-up
            tech_nodes, zone_nodes = set(), set()# remember: possibly multiple
            for zone in zones:
                hvac = zone.inputs(MODE)# MODE = skt.identifier :-)
                for link in hvac.links:
                    tech_nodes.add(link.from_node); break
                else: zone_nodes.add(zone)
            ### check for tech-nodes, which are linked to other zone-nodes too.
            #   we don't want to accidentally delete or override them; instead
            #   we simply unlink our zones from them and update our references.
            # NOTE: we only remove/add when necessary (invasiveness).
            # NOTE: store links first, remove them later all-at-once.
            # FIXME: how to deal with the vi-suite's update-bogging ?
            for tnd in tech_nodes:
                links, all_links = set(), set(tnd.links)
                for lnk in all_links: # just checking links
                    if lnk.to_node in zones: links.add(lnk)
                if  links!=all_links: # linked to other zone-nodes => Unlink us!
                    for lnk in links: # if links was empty, tnd hadn't been here
                        zone_nodes.add(lnk.to_node)
                        envi_network.links.remove( lnk )
                    tech_nodes.remove(tnd)  # this node was successfuly unlinked
            ### now unlink/link up, and add our own node if necessary:
            TND = None # dummy, in case we are unlinked (None=>dont update info)
            if add:
                ## fetch our own linked technode ...
                for skt in node.inputs[idx_a:idx_b]:# going through the sockets
                    if skt.bl_idname==skt_bl_idname:# there's at most one
                        for lnk in skt.links: TND = lnk.from_node; break
                ## add target technode if needed ...
                if not tech_nodes:
                    tnd = envi_network.nodes.new(tnd_bl_idname)
                    tnd.location = X0-400, Y0-idx_m*60# stacked left of the zones
                    tnd.hide = True
                    tech_nodes.add(tnd)
            ## link and update, or remove the technode
            for i, tnd in enumerate(tech_nodes):
                if add: # we mean to add/update
                    if i==0:# just link all there is to link to the first one...
                        for zone in zone_nodes:
                            envi_network.links.new(# MODE = skt.identifier :-)
                                tnd.outputs[MODE], zone.inputs[MODE])
                    # supply with data from our own linked technode (just copy
                    if TND: _copy_settings(TND, tnd, settings)# our settings)
                # not adding? we want to remove then
                else: envi_network.nodes.remove(tnd)
            # ## ### Finished, i think?


class NODETREE_OT_goto_group(bpy.types.Operator):
    'Enter this node group. If append=False set tree to JUST THIS node group.'
    bl_idname = 'nodetree.goto_group'
    bl_label = 'Go To Group'

    tree: bpy.props.StringProperty(default="")
    append: bpy.props.BoolProperty(default=True)# (basically toggle into group)

    def execute(self, context):
        tree = bpy.data.node_groups.get(self.tree, None)
        if not tree:self.report({'WARNING'},"Tree not found: "+self.tree)
        elif self.append: context.space_data.path.append(tree)#push group=>stack
        else:# set tree as active tree and space_data.path to only this tree
            context.space_data.tree_type = tree.bl_idname#ensure tree_type first
            context.space_data.node_tree = tree# space_data.path=tree & activate
        return {('FINISHED','ERROR')[not tree]}

# TODO: unfinished!!!
class IFC2SIM_OT_setup_tree(bpy.types.Operator):
    'Set up all nodes needed for the IFC2SIM-Workflow'
    bl_idname = 'ifc2sim.setup_tree'
    bl_label = 'Set up tree'
    bl_options= {"REGISTER"}#,"UNDO"}# we might be too large for the undo-stack?
    bl_description = "Setup & Run IFC2SIM project."

    filepath: bpy.props.StringProperty(subtype="FILE_PATH")
    filter_glob: bpy.props.StringProperty(default="*.ifc;*.ifczip;*.ifcxml", options={"HIDDEN"})
    is_advanced: bpy.props.BoolProperty(name="Enable Advanced Mode", default=False)
    ## skip any actual actions. Only add and link the nodes
    only_add: bpy.props.BoolProperty("only add", default=False)# not implemented

    # NOTE: we pull these from IFC2SIM_OT_import_ifc so we can start the tree
    #       because the import-node only has an output to work with after import
    invoke = IFC2SIM_OT_import_ifc.invoke#=> allow for importer-interface to run
    #import = IFC2SIM_OT_import_ifc.execute#=> the actual import

    #def invoke(self, context, event):## we dont need to implement 'invoke' here
        #context.window_manager.fileselect_add(self)## we can just steal it
        #self.node = context.node if hasattr(context, "node") else None
        #return {"RUNNING_MODAL"}

    def execute(self, context):
        self.only_add = True # we're not yet finished ... WARNING: UNSTABLE CODE
        return self.execute_real(context)
    # real
    def execute_real(self, context):
        tree    = context.space_data.edit_tree
        nodes   = tree.nodes
        offx, offy = 200, 200 # base-offset for our nodes, to set their position
        ### These are the nodes we need. We will only instance them once needed.
        #   NodeImportIFC, 2*NodeObjectsFilter, NodeSetMaterials,
        #   NodeCutOpenings, NodeExportEnViGeometry and NodeTechAdd
        # NOTE: We instance just-in-time to ease garbage-collection
        ### add nodes, link them up and apply operators
        #
        #=# first let's do the import stuff: NOTE: this we needed the invoke for
        imp_ifc                 = nodes.new("NodeImportIFC")
        all_objs = imp_ifc.outputs["All Objects"]# we have exactly one output.
        if not self.only_add: ### if real
            VL = context.view_layer; objects = list(VL.objects)# back-up the objects
            fn = bpy.ops.bim.load_project('EXEC_DEFAULT', filepath=self.filepath,
                filter_glob=self.filter_glob, is_advanced=self.is_advanced)
            # store the new objects in an appropriate output
            #if "All Objects" in imp_ifc.outputs:
            #    all_objs = imp_ifc.outputs["All Objects"]
            #else: all_objs=imp_ifc.outputs.new('SktObjects', "All Objects")
            all_objs.objects = [obj for obj in VL.objects if obj not in objects]
        #
        #=# next extract zones, windows and doors:
        split_objs              = nodes.new("NodeObjectsFilter")
        split_objs.location.x   =    1*offx
        tree.links.new(all_objs, split_objs.inputs[0])# link the input
        # add three filters and their outputs (doesn't actually need a context):
        # NOTE: part could be done with .upd_flt_add(None), which seems overkill
        for name in ["Zones","Doors","Windows"]:
            split_objs.outputs.new('SktObjects',name)   # socket skt
            flt = split_objs.filters_collection.add()   # filter flt
            flt.name = name
            flt.filter_type="flt_IFCType"
            if   name=="Doors"  : flt.flt_IFCType="Opening/Door"
            elif name=="Windows": flt.flt_IFCType="Opening/Window"
            else                : flt.flt_IFCType="IfcSpace" # handle zones
        split_objs.index = 2
        if not self.only_add: ### if real
            # apply filters, storing the objects in their respective outputs:
            split_objs.upd_update_all(None)# doesn't actually need any context
        #
        #=# next apply some general materials (assume&create them if necessary):
        set_mat     = nodes.new("NodeSetMaterials")
        set_mat.location        =    2*offx, 1*offy
        tree.links.new(split_objs.outputs[0], set_mat.inputs[0])# link the input
        if not self.only_add: ### if real
            def fetch_materials(name):
                if name in bpy.data.materials: return bpy.data.materials[name]
                return bpy.data.materials.new(name)
            set_mat.material_ceiling    = fetch_materials("Ceiling")
            set_mat.material_floor      = fetch_materials("Floor")
            set_mat.material_north      = fetch_materials("North")
            set_mat.material_south      = fetch_materials("South")
            set_mat.material_east       = fetch_materials("East")
            set_mat.material_west       = fetch_materials("West")
            set_mat.upd_set_materials(None)
        #
        #=# after the general materials, add some doors and windows
        cut_openings= nodes.new("NodeCutOpenings")
        cut_openings.location   =    2*offx,  -0.25*offy
        tree.links.new(split_objs.outputs[0], cut_openings.inputs[0])#link input
        tree.links.new(split_objs.outputs[1], cut_openings.inputs[1])#link input
        tree.links.new(split_objs.outputs[2], cut_openings.inputs[2])#link input
        if not self.only_add: ### if real
            cut_openings.material_door   = mat_door   = fetch_materials("Door")
            cut_openings.material_window = mat_window = fetch_materials("Window")
            #fn = bpy.ops.ifc2sim.cut_openings('EXEC_DEFAULT')# cut in our openings
            # try with direct context override
            #fn = bpy.ops.ifc2sim.cut_openings({"node":cut_openings})
            # made it a core function?
            # inlined it!
        zones   = split_objs.outputs[0].objects
        doors   = split_objs.outputs[1].objects
        windows = split_objs.outputs[2].objects
        if not self.only_add: ### if real
            if zones:
                collection = zones[0].users_collection[0]
                if doors:
                    #print("CUTTING DOORS ####################################")
                    zones, cutter = cut_blocks_new(zones, doors,   collection,
                                                cut_openings.material_door)
                    bpy.data.objects.remove(cutter)# cleanup-job: delete cutters
                if windows:
                    #print("CUTTING WINDOWS ##################################")
                    zones, cutter = cut_blocks_new(zones, windows, collection,
                                                cut_openings.material_window)
                    bpy.data.objects.remove(cutter)# cleanup-job: delete cutters
        skt = cut_openings.outputs.get("Zones")# get output socket, if available
        if not skt: skt=cut_openings.outputs.new('SktObjects', "Zones")#else add
        if not self.only_add: ### if real ## muss das echt bei dummies weg?
            skt.objects = zones                    # push the modified zones.
        #
        #=# next set up the zones to be zones and then export them
        exp_geo                 = nodes.new("NodeEnViGeometryExport")# NodeEnViGeometryExport # NodeExportEnViGeometry
        exp_geo.location        =    3*offx, 1*offy
        tree.links.new(cut_openings.outputs[0], exp_geo.inputs[0])#link up input
        if not self.only_add: ### if real
            exp_geo.upd_setup_zones(None)# zones-setup doesn't actually need context
            # for whatever reason, engexport uses invoke instead of execute?
            fn = bpy.ops.node.engexport({"node": exp_geo}, 'INVOKE_DEFAULT')
            #
            #=# link up the exported zones and fix their volumes. NOTE: The vi-suite
            link_exported_zones()   # inherits this bug from bmesh, when calculating
            set_correct_zone_volume()# the volume of non-triangulated bodies with it.
        #
        #=# next extract the zones into their respective usecases
        split_objs2             = nodes.new("NodeObjectsFilter")
        split_objs2.location    =    3*offx, 60
        tree.links.new(cut_openings.outputs[0], split_objs2.inputs[0])# link up
        # add three filters and their outputs (doesn't actually need a context):
        # NOTE: part could be done with .upd_flt_add(None), which seems overkill
        for name in ["Bath","Hall","Living","Sleeping"]:
            split_objs2.outputs.new('SktObjects',name)  # socket skt
            flt = split_objs2.filters_collection.add()  # filter flt
            flt.name = name
            flt.filter_type="flt_ZoneType"
            flt.flt_ZoneType=name.upper()#"BATH" HALL LIVING SLEEPING
        split_objs2.index = 3
        if not self.only_add: ### if real
            # apply filters, storing the objects in their respective outputs:
            split_objs2.upd_update_all(None)# doesn't actually need any context
        #
        #=# now we're ready to set-up our technodes, (carefull, it may crash us)
        add_technode            = nodes.new("NodeTechAdd")
        add_technode.location.x =    5*offx
        # next add four zones-inputs and link them up
        # NOTE: part could be done with .upd_ptr_add(None), which seems overkill
        for name in ["Bath", "Hall", "Living", "Sleeping"]:
            skt = add_technode.inputs.new('SktObjects', name)   # socket    skt
            tnd = add_technode.technodes.add()                  # technode  tnd
            tnd.name = name
            skt.show_expanded=False
        tree.links.new(split_objs2.outputs[0], add_technode.inputs[0])# link up
        tree.links.new(split_objs2.outputs[1], add_technode.inputs[1])# link up
        tree.links.new(split_objs2.outputs[2], add_technode.inputs[2])# link up
        tree.links.new(split_objs2.outputs[3], add_technode.inputs[3])# link up
        if not self.only_add: ### if real
            # now we need sensible dummy data. This is where we break for now.
            # TODO: continue setup
            #{bath 24 deg, hall 18deg, Sleeping 18, Living 21}
            ###
            pass
        return {"FINISHED"}

# TODO: Object.set.material vs faces from material vs faces.set.material
###


# END ##########################################################################
###                           !!! REGISTRATION !!!                           ###
#######################################80############################### BEGIN #
### Node Categories ###
# Node categories are a python system for automatically
# extending the Add menu, toolbar panels and search operator.
# For more examples see release/scripts/startup/nodeitems_builtins.py

import nodeitems_utils
from nodeitems_utils import NodeCategory, NodeItem

# our own base class with an appropriate poll function,
# so the categories only show up in our own tree type
class MyNodeCategory(NodeCategory):
    @classmethod
    def poll(cls, context):
        return context.space_data.tree_type == IFC_SIM_Tree.bl_idname

# Generic node-operator
class NodeItemOperator:
    def __init__(self, op_id, poll=None, **settings):
        # settings: (text, text_ctxt, translate, icon, emboss, depress, icon_value)
        self.draw=lambda s,col,ctxt:col.operator(op_id, **settings)
        self.poll=poll


### Generic sub-category. Can be nested. Some bug with the poll-function prevents
#   use in the top-layer (in place of NodeCategory). FIXME: fix polling-issues
#   NOTE: One might have to subclass from this to truely add a poll-function?!?
class SubCategory:
    def __init__(self, identifier, name, description=None, poll=None, items=None):
        ### first check wether we already exist:
        idname = "NODE_MT_subcategory_" + identifier
        if idname in dir(bpy.types): pass
            #raise KeyError("Node subcategory '%s' already registered" % identifier)
        ### No? Good!
        # Set up some variables:
        self.identifier, self.name, self.description = identifier, name, description
        _poll = poll if callable(poll) else lambda cls,context:True
        self.poll = None#classmethod(lambda cls,context:_poll(self,context)) # no sense in polling here
        # works as draw function for menus # stolen from nodeitems_utils without
        def draw_node_item(self, context): # any modifications.
            layout = self.layout
            col = layout.column()
            for item in self.category.items(context):
                item.draw(item, col, context)
        # add and register ourselves a menu type
        self.menu_type = type(idname, (bpy.types.Menu,), {
            "bl_space_type": 'NODE_EDITOR',
            "bl_label": name,
            "category": self,          # we are the category
            "poll": classmethod(_poll),# we define poll ourselves
            "draw": draw_node_item,
        })
        bpy.utils.register_class(self.menu_type)
        # Set up self.items # stolen from nodeitems_utils without modifications.
        if items is None:
            self.items = lambda context: []
        elif callable(items):
            self.items = items
        else:
            def items_gen(context):
                for item in items:
                    if item.poll is None or context is None or item.poll(context):
                        yield item
            self.items = items_gen
    @staticmethod
    def draw(self, layout, context):
        # just like in nodeitems_utils, we draw the menu, but there's no sense in polling
        layout.menu("NODE_MT_subcategory_%s" % self.identifier)
    #@classmethod
    #def _poll(cls,context): return False
    #@classmethod
    #def poll(cls, context): return True#cls._poll(cls,context)

classes = (
    ### UI-Lists
    POINTERS_UL_general,
    FILTERS_UL_simple,
    TECHNODE_UL_simple,

    IFC_SIM_Tree,

    MyCustomSocket,
    SktObjects,
    SktCollections,

    SocketNode,
    SocketTesterNode,
    SocketValueTesterNode,
    ### Nodes
    #   for objects
    NodeObjectsFilterLive,
    NodeObjectsJoin,
    NodeObjectsInput,
    NodeObjectsViewer,
    #   for collections
    NodeCollectionsInput,
    NodeCutOpenings,
    NodeSetMaterials,
    NodeEnViGeometryExport, ## NodeExportEnViGeometry, NodeLinkExportedZones, NodeSetVolume,
    NodeImportIFC,
    NodeObjectsFilter,
    NodeCollectionsAddObjects,
    NodeTechAdd,
    NodeAnalyseZones,
    NodeEnviContext,

    ### Operators
    IFC2SIM_OT_cut_openings,
    IFC2SIM_OT_import_ifc,
    IFC2SIM_OT_collection_items,
    IFC2SIM_OT_update_technodes,
    NODETREE_OT_goto_group,# universal: toggle into / enter a certain node group
    IFC2SIM_OT_setup_tree,
)

### currently not used. This adds nodes as NodeItems and appends them to classes if necessary.
def NodeItemList(*nodes):
    classes[:]=dict.fromkeys(classes+nodes)# since PY3.7 dicts are order-preserving
    return [NodeItem(node.bl_idname) for node in nodes]

# all categories in a list
node_categories = [
    # identifier, label, items list
    MyNodeCategory('TESTNODES', "Some Nodes for testing", items=[
        # our basic node
        NodeItem("SocketTesterNode"),
        NodeItem("SocketNode"),
        NodeItem("SocketValueTesterNode"),
        NodeItem("NodeCutOpenings"),
        # prepare, export and link zones, set their volume
        NodeItem("NodeEnViGeometryExport"), ## NodeItem("NodeExportEnViGeometry"), NodeItem("NodeLinkExportedZones"), NodeItem("NodeSetVolume"),
        NodeItem("NodeImportIFC"),
        NodeItem("NodeTechAdd"),
        NodeItem("NodeAnalyseZones"),
        NodeItem("NodeEnviContext"),
        NodeItemOperator("IFC2SIM_OT_setup_tree"),
    ]),
    MyNodeCategory('OBJECTS', "Object", items=[
        #NodeItem("NodeObjectsFilterLive"),
        NodeItem("NodeObjectsFilter"),
        NodeItem("NodeObjectsJoin"),
        NodeItem("NodeObjectsInput"),
        NodeItem("NodeObjectsViewer"),
    ]),
    MyNodeCategory('MATERIAL', "Material", items=[
        NodeItem("NodeSetMaterials"),
    ]),
    MyNodeCategory('COLLECTION', "Collection", items=[
        NodeItem("NodeCollectionsInput"),
        NodeItem("NodeCollectionsAddObjects"),
    ]),
    MyNodeCategory('LAYOUT', "Layout", items=[
        NodeItem("NodeFrame"),
        NodeItem("NodeReroute")
    ]),
    ### Nestable pseudo-NodeCategory. Needs real NodeCategory above it.
    #    SubCategory('TestSubCatAB', 'test subcat', items=[])
    #MyNodeCategory('VIEWERS', "Viewers", items=[NodeItem(idname) for idname in (
        #"NodeObjectsFilterLive","NodeObjectsJoin","")]),
    #MyNodeCategory('INPUTS',  "Inputs",  items=[NodeItem(idname) for idname in (
        #"NodeObjectsInput", "NodeCollectionsInput")]),
    #MyNodeCategory('OTHERNODES', "Other Nodes", items=[
    #    # NodeItems can have additional 'settings', to apply them to new nodes.
    #    # NB: settings values are stored as string expressions,
    #    # for this reason they should be converted to strings using repr()
    #    NodeItem("SocketTesterNode", label="Node A", settings={
    #        "my_string_prop": repr("Lorem ipsum dolor sit amet"),
    #        "my_float_prop": repr(1.0)}),
    #]),
]


def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)

    nodeitems_utils.register_node_categories('CUSTOM_NODES', node_categories)


def unregister():
    nodeitems_utils.unregister_node_categories('CUSTOM_NODES')

    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)


if __name__ == "__main__":
    try:
        unregister()
    except:
        pass
    register()

# END ##########################################################################
###                              !!! NOTES !!!                               ###
#######################################80############################### BEGIN #
#   Fuer Aufbereitung-der-Zonen benoetige ich einen Filter fuer IfcOpenings    #
#       und einen Filter fuer Space simple (unser Typ fuer Zonen)              #
#   Import-Node funktioniert erst mal.  Filter-Node könnte erweitert werden.   #
# + Vielleicht einen Object-Split-By node einbauen? Fertig (vorerst)!          #
#   NOTE: Object-Split-By Node unterstützt Filter für Namen, IFC-Typ, Zonentyp,#
#       und Position (Origin, BoundingBox, Vertex) weitgehend oder teilweise.  #
#       Die Darstellung der Outputs in der UI-List macht name_dummy überflüssig#
#   NOTE: The socket has an identifier (which is set fix on creation), but     #
#       the search within inputs/outputs goes by name. If we create a string   #
#       (single-use) property we may search for the identifier ourselves.      #
#   TODO: Compare polygons for belonging to the same wall-segment:             #
#       A) compare angle B) Hesse Normal form => distance (distance of centers?#
#       C) project vertices (onto the mean plane?) D) 2D-intersect             #
#   Wären mengen-operationen interessant? Concatenate/Union/Difference/etc.    #
#   NOTE: der Chart-Node erzeugt für die x und jede der drei Y-achsen einen    #
#       eigenen Socket. der unterscheidet sich aber eigentlich NUR darin,      #
#       welche items von dem jeweiligen Enum verwendet werden => Wie machen wir#
#       das bei uns besser? => Dem Chart einfach nur einen Socket geben und die#
#       items als items-funktion übergeben. DAFÜR GIBT ES DIE!!!               #
#       ODER noch besser: Wir machen einen einzigen Socket, der beliebig viele #
#       y-achsen gestattet. => Update-Function und eine collection an enums, zB#
#   NOTE: Im "IFC2SIM_OT_setup_tree" wurden die execute-funktion mit und ohne  #
#       reale Aktionen zusammengeführt. Aktuell war ohnehin nur die Version    #
#       ohne reale Aktion verfügbar. Das ist jetzt auch so.                    #
#                                                                              #
################################################################################
# END ##########################################################################
