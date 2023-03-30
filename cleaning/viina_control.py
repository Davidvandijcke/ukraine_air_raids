from arcgis import GIS
from arcgis.features import FeatureSet
from arcgis.mapping import WebMap

gis = GIS()

the_map = WebMap(gis.content.get('https://experience.arcgis.com/experience/6314230485b64d71b3b9ef65192a205c'))

for l in the_map.layers[:-1]:
    fs = FeatureSet.from_dict(l['featureCollection']['layers'][0]['featureSet'])