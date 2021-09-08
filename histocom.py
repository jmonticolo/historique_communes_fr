#!/usr/bin/env python3
import os
import uuid
from collections import Counter
from datetime import date, datetime, timedelta
from math import log
from os.path import isdir, isfile, join

import geopandas as gam
import pandas as am
import pytz
from dateutil.relativedelta import relativedelta
from pyproj import Transformer
from shapely.geometry import MultiPolygon, Point, Polygon
from shapely.ops import unary_union


class GeomObj(object):
    """Geometric object class"""

    def __init__(self):
        self._gid = str(uuid.uuid4())
        self._attrs = {}
        self._shape = None

    @property
    def gid(self):
        return self._gid

    @property
    def attributes(self):
        return self._attrs

    def attribute(self, attribut: str):
        if attribut in self._attrs:
            return self._attrs[attribut]
        else:
            return None

    def set_attribute(self, attribut: str, value):
        self._attrs[attribut] = value

    def copy_attributes(self, other_geo, attributes: list):
        self._attrs.update({attr: other_geo.attribute(attr) for attr in attributes})

    @property
    def shape(self):
        return self._shape

    @shape.setter
    def shape(self, geometry):
        self._shape = geometry

    def has_geom(self):
        return self._shape is not None

    def clean_geom(self):
        """Nettoie la geometrie"""
        # si la geometrie est deja nettoyee
        if self.attribute("is_clean"):
            return self

        buf = 0.000001
        # enleve les artefacts tres fins
        self._shape = self._shape.buffer(-buf).buffer(buf)
        # enleve les trous interieurs
        if self._shape.geom_type == "Polygon":
            self._shape = unary_union(Polygon(self._shape.exterior))
        elif self._shape.geom_type == "MultiPolygon":
            parts = []
            for poly in self._shape:
                parts.append(Polygon(poly.exterior))

            self._shape = unary_union(parts)

        self.set_attribute("is_clean", True)


class ComObj(object):

    P_START, P_END = range(2)

    def __init__(
        self,
        admin_code: str = None,
        mod: int = None,
        name: str = None,
        period_start=None,
        period_end=None,
        prev_coms: set = set(),
        next_coms: set = set(),
        population: int = None,
        geom: GeomObj = None,
    ):
        self._admin_code = admin_code
        self._mod = mod
        self._name = self.correct_name(name)
        self._period_start = period_start
        self._period_end = period_end
        self._prev_coms = prev_coms.copy()
        self._next_coms = next_coms.copy()
        self._pop = population
        if geom is None:
            self._geom = GeomObj()
        else:
            self._geom = geom

        if self._geom.attribute("source") is None:
            self._geom.set_attribute("source", "INSEE")

    @property
    def admin_code(self) -> str:
        return self._admin_code

    @property
    def mod(self) -> int:
        return self._mod

    @mod.setter
    def mod(self, value: int):
        self._mod = value

    @staticmethod
    def correct_name(text: str) -> str:
        correct = {"¼": "Œ", "½": "œ"}
        for old, new in correct.items():
            text = text.replace(old, new)

        return text

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, value: str):
        self._name = self.correct_name(value)

    @property
    def period_start(self):
        return self._period_start

    @period_start.setter
    def period_start(self, value):
        self._period_start = value

    @property
    def period_end(self):
        return self._period_end

    @period_end.setter
    def period_end(self, value):
        self._period_end = value

    @property
    def prev_coms(self) -> set:
        return self._prev_coms

    @prev_coms.setter
    def prev_coms(self, value: set):
        self._prev_coms = value.copy()

    def prev_coms_add(self, value: str):
        self._prev_coms.add(value)

    @property
    def next_coms(self) -> set:
        return self._next_coms

    @next_coms.setter
    def next_coms(self, value: set):
        self._next_coms = value.copy()

    def next_coms_add(self, value: str):
        self._next_coms.add(value)

    @property
    def population(self) -> int:
        return self._pop

    @population.setter
    def population(self, value):
        if isinstance(value, float):
            self._pop = value * 1000
        elif isinstance(value, int):
            self._pop = value
        else:
            print(type(value))

    @property
    def geom(self) -> GeomObj:
        return self._geom

    @geom.setter
    def geom(self, value: GeomObj):
        self._geom = value


class Histocom(object):
    """Build the french cities history with INSEE and IGN data"""

    def __init__(self):
        self.communes = {}
        # code du champ mod: priorite de traitement
        self._mods = {
            10: 1,  # Changement de nom
            20: 2,  # Création
            21: 2,  # Rétablissement
            30: 3,  # Suppression
            31: 5,  # Fusion simple
            32: 7,  # Création de commune nouvelle
            33: 6,  # Fusion association
            34: 1,  # Transformation de fusion association en fusion simple
            41: 1,  # Changement de code dû à un changement de département
            50: 1,  # Changement de code dû à un transfert de chef-lieu
            70: 9,  # Transformation de commune associé en commune déléguée
        }
        self.insee_df = None
        self.date_start = date(1943, 1, 1)
        self.insee_period = None
        self._histo_df = None
        self.mergers = {}
        self.splits = {}
        self.timezones = {
            "paris": pytz.timezone("Europe/Paris"),
            "971": pytz.timezone("America/Guadeloupe"),
            "972": pytz.timezone("America/Martinique"),
            "973": pytz.timezone("Indian/Mayotte"),
            "974": pytz.timezone("Indian/Reunion"),
            "976": pytz.timezone("America/Cayenne"),
            "utc": pytz.utc,
        }

    @property
    def mods(self):
        return sorted(self._mods.keys())

    def get_mod(self, insee_code, period, type_period):
        """
        Permet l'acces a une modification pour une commune
        en specifiant une date et le type de date
        """
        if insee_code not in self.communes:
            return None

        if type_period not in [ComObj.P_START, ComObj.P_END]:
            return None

        for m in self.communes[insee_code]:
            if type_period == ComObj.P_START:
                if m.period_start == period:
                    return m
            elif type_period == ComObj.P_END:
                if m.period_end == period:
                    return m

        return None

    def priority(self, mod: int):
        return self._mods[mod]

    def build_insee(self, path: str):
        """Construit l'historique des communes a partir du fichier INSEE"""
        df = am.read_csv(path)
        cols = ["mod", "date_eff"]
        for i in ["typecom", "com", "libelle"]:
            for j in ["_av", "_ap"]:
                cols.append(i + j)

        if not all([i in df for i in cols]):
            raise ValueError("Missing column")

        file_mods = sorted(df["mod"].unique())
        if not file_mods == self.mods:
            missing_mods = [m for m in self.mods if m not in file_mods]
            raise ValueError(f"Mod(s) {', '.join(missing_mods)} not taken into account")

        # ajoute la serie priorite par mod
        df["mod_priority"] = df.apply(lambda row: self.priority(row["mod"]), axis=1)
        # convertit la date
        df.date_eff = df.apply(
            lambda row: str(
                datetime.strptime(row["date_eff"], "%d/%m/%y").date()
                - relativedelta(years=100)
                if datetime.strptime(row["date_eff"], "%d/%m/%y").date()
                > datetime.now().date()
                else datetime.strptime(row["date_eff"], "%d/%m/%y").date()
            ),
            axis=1,
        )
        df = df.sort_values(by=["date_eff", "mod_priority", "com_ap", "com_av"])

        self.insee_df = df.loc[
            lambda df: (df["typecom_av"] == "COM") & (df["typecom_ap"] == "COM")
        ]

        dates_eff = self.insee_df["date_eff"].unique()
        for d in dates_eff:
            self.date_modif(self.insee_df.loc[lambda df: df["date_eff"] == d])

        self.insee_period = [self.date_start, date(*[int(v) for v in d.split("-")])]

    def date_modif(self, dataframe):
        """
        dataframe type: Pandas dataframe
        """
        modifs = {20: {}, 21: {}, "3_0123": {}}
        day = timedelta(days=1)
        for data in dataframe.itertuples():
            date_eff = date(*[int(v) for v in data.date_eff.split("-")])
            mod = "3_0123" if data.mod in [30, 31, 32, 33] else data.mod
            com_av = data.com_av
            com_ap = data.com_ap
            lib_av = data.libelle_av
            lib_ap = data.libelle_ap
            com_pop = False

            # inutiles
            if mod == 34 and com_av == com_ap and lib_av == lib_ap:
                continue
            elif mod == 70:
                continue

            # ajout de l'ancienne commune si non presente
            if com_av not in self.communes:
                self.communes[com_av] = [
                    ComObj(
                        admin_code=com_av,
                        mod=data.mod,
                        name=lib_av,
                        period_start=self.date_start,
                        period_end=date_eff - day,
                        next_coms={com_ap},
                    )
                ]

            prev_com = self.communes[com_av][-1]
            next_com = ComObj(
                admin_code=com_av,
                mod=data.mod,
                name=lib_ap,
                period_start=date_eff,
                prev_coms={com_av},
            )

            if mod in [10, 34, 41, 50]:
                # ajout de la nouvelle commune
                next_com.geom = prev_com.geom
                if com_ap in self.communes:
                    self.communes[com_ap].append(next_com)
                else:
                    self.communes[com_ap] = [next_com]
            elif mod in [20, 21, "3_0123"]:
                if com_ap not in modifs[mod]:
                    modifs[mod][com_ap] = next_com
                else:
                    modifs[mod][com_ap].prev_coms_add(com_av)

                # preparation des operations geometriques
                if data.mod == 21:
                    skey = prev_com.geom
                    sval = modifs[mod][com_ap].geom
                    if skey in self.splits:
                        self.splits[skey].append(sval)
                    else:
                        self.splits[skey] = [sval]
                elif data.mod in [31, 32, 33]:
                    # print(modifs)
                    mkey = modifs[mod][com_ap].geom
                    mval = prev_com.geom
                    if mkey in self.mergers:
                        self.mergers[mkey].append(mval)
                    else:
                        self.mergers[mkey] = [mval]

                if mod == "3_0123":
                    # modifications le meme jour
                    for m in [20, 21]:
                        if com_av in modifs[m]:
                            modifs[m].pop(com_av)
                            com_pop = True

            # mise a jour des anciennes communes
            if prev_com.period_end is None:
                prev_com.period_end = date_eff - day
            if not com_pop:
                prev_com.next_coms_add(com_ap)
            if (
                mod == "3_0123"
                and com_ap in self.communes
                and self.communes[com_ap][-1].period_end is None
            ):
                self.communes[com_ap][-1].period_end = date_eff - day

            # modifications le meme jour
            if prev_com.period_start == date_eff:
                df = dataframe.loc[
                    lambda df: (df["date_eff"] == str(date_eff))
                    & (df["com_ap"] == com_av)
                    & (df["libelle_ap"] == lib_av)
                ]
                self.communes[df["com_av"].values[0]][-1].next_coms_add(com_ap)
                self.communes[com_av].pop(-1)
                if len(self.communes[com_av]) == 0:
                    self.communes.pop(com_av)

        # ajout des nouvelles communes
        for m, d in modifs.items():
            for k, v in d.items():
                if k in self.communes:
                    self.communes[k].append(v)
                else:
                    self.communes[k] = [v]

    @staticmethod
    def get_date(date_str: str):
        date_parts = date_str.split("/")[-1].split("-")
        year = int(date_parts[0])
        month = int(date_parts[1]) if len(date_parts) == 2 else 1
        return date(year, month, 1)

    def insert_month_dir(self, path_list: list, new_path) -> list:
        new_date = self.get_date(new_path)
        for i, d in enumerate(path_list):
            if new_date > self.get_date(d):
                path_list.insert(i, new_path)
                break

        return path_list

    def build_geom(self, path: str):
        """
        Mets a jour l'historique au niveau geometrique et complete
        avec les communes n'ayant jamais eu de modifications
        """
        if isdir(path):
            geom_dirs = [f for f in os.scandir(path) if isdir(f)]
            # annees
            paths = sorted(
                [d.path for d in geom_dirs if len(d.name) == 4], reverse=True
            )
            # mois
            geo_month = sorted(
                [d.path for d in geom_dirs if len(d.name) == 7], reverse=True
            )
            # ajout des mois aux annees
            for p in geo_month:
                paths = self.insert_month_dir(paths, p)

            self._histo_df = self.df()
            self._histo_df["has_geom"] = None
            for geo_path in paths:
                print(geo_path)
                self.geom_insee(geo_path)

            maj = True
            while maj:
                maj = False

                m1 = True
                while m1:
                    m1 = False
                    m1 = self.update_same_geom()
                    maj = True if m1 else maj

                m1 = True
                while m1:
                    m1 = False
                    m2 = self.merge_geom(self.mergers)
                    m3 = self.merge_geom(self.splits)
                    m1 = any([m2, m3])
                    maj = True if m1 else maj

                m1 = True
                while m2:
                    m1 = False
                    m2 = self.geom_by_difference(self.mergers)
                    m3 = self.geom_by_difference(self.splits)
                    m1 = any([m2, m3])
                    maj = True if m1 else maj

            self.geom_after_insee(paths)
            self.update_latlon()

    @staticmethod
    def read_source(path: str):
        if isfile(join(path, "source.txt")):
            with open(join(path, "source.txt"), "r") as f:
                source_file = f.readline().strip()
                if source_file[:13] == "ADMIN-EXPRESS":
                    return "AdminExpress"
                elif source_file[:6] == "GEOFLA":
                    return "Geofla"
                else:
                    return None
        else:
            return None

    def geom_insee(self, path: str):
        """
        Recupere les geometries et attributs pour les communes presentes
        dans le fichier INSEE
        """
        dirdate = self.get_date(path)
        # hors periode INSEE
        if dirdate > self.insee_period[1]:
            return

        geodf = gam.read_file(join(path, "COMMUNE.shp"))
        source = self.read_source(path)

        chfdf = None
        if isfile(join(path, "CHEF_LIEU.shp")):
            chfdf = gam.read_file(join(path, "CHEF_LIEU.shp"))

        df2 = self._histo_df.loc[
            lambda df: df["has_geom"].isnull()
            & (df["period_start"] <= dirdate)
            & ((dirdate <= df["period_end"]) | df["period_end"].isnull())
        ]

        for t in df2.itertuples():
            g2 = geodf.loc[lambda df: df["INSEE_COM"] == t.insee_code]
            if g2["INSEE_COM"].count() == 1:
                record = g2.iloc[0]
                if not t.geom.has_geom():
                    if "ID_GEOFLA" in record:
                        t.geom.set_attribute("id", record.ID_GEOFLA)
                    elif "ID" in record:
                        t.geom.set_attribute("id", record.ID)

                    if "POPULATION" in record:
                        com_mod = self.get_mod(
                            t.insee_code, t.period_start, ComObj.P_START
                        )
                        pop = record.POPULATION
                        # jusqu'en 2013, la population est en milliers d'individus
                        if pop.dtype.name == "float64":
                            pop = pop * 1000

                        com_mod.population = int(pop)

                    x_wgs, y_wgs = None, None
                    if "X_CHF_LIEU" in record:
                        project = Transformer.from_crs("EPSG:2154", "EPSG:4326")
                        x_l93, y_l93 = record.X_CHF_LIEU, record.Y_CHF_LIEU
                        if y_l93 < 1000000:
                            factor = 6 - int(log(y_l93, 10))
                            x_l93 = x_l93 * 10 ** factor
                            y_l93 = y_l93 * 10 ** factor

                        y_wgs, x_wgs = project.transform(x_l93, y_l93)
                    elif chfdf is not None:
                        c2 = chfdf.loc[lambda df: df["INSEE_COM"] == t.insee_code]
                        if c2["INSEE_COM"].count() == 1:
                            x_wgs, y_wgs = c2.iloc[0].geometry.x, c2.iloc[0].geometry.y

                    t.geom.set_attribute("latitude", y_wgs)
                    t.geom.set_attribute("longitude", x_wgs)
                    t.geom.set_attribute("source", source)
                    t.geom.shape = record.geometry
                    t.geom.set_attribute("is_clean", False)
                    self._histo_df.loc[t.Index, "has_geom"] = True
            elif g2["INSEE_COM"].count() > 1:
                raise ValueError(f"Duplicate on INSEE_COM found: {t.insee_code}")

    def update_same_geom(self) -> bool:
        """Mets a jour la geometrie pour des communes ayant la meme geometrie"""
        up = False
        for k, v in self.communes.items():
            for i in range(len(v) - 1, -1, -1):
                if v[i].mod in [31, 32, 33]:
                    # s'il n'y a eu qu'une seule modification depuis la fusion
                    if len(v[i:]) == 2:
                        if (
                            v[1].mod == 21
                            and len(v[0].next_coms) == 1
                            and v[0].next_coms == v[1].prev_coms
                        ):
                            if v[0].geom.has_geom() and not v[1].geom.has_geom():
                                v[1].geom.shape = v[0].geom.shape
                                v[1].geom.copy_attributes(
                                    v[0].geom,
                                    ["is_clean", "latitude", "longitude", "source"],
                                )
                                up = True
                            elif not v[0].geom.has_geom() and v[1].geom.has_geom():
                                v[0].geom.shape = v[1].geom.shape
                                v[0].geom.copy_attributes(
                                    v[1].geom,
                                    ["is_clean", "latitude", "longitude", "source"],
                                )
                                up = True
                    else:
                        if v[i].next_coms == {k}:
                            lc = [k]
                            nex = set()
                            for m in v[i + 1 :]:
                                if m.mod == 21 and m.prev_coms == {k}:
                                    lc = [n for n in lc if n not in nex]
                                    if lc == []:
                                        if (
                                            v[i].geom.has_geom()
                                            and not m.geom.has_geom()
                                        ):
                                            m.geom.shape = v[i].geom.shape
                                            m.geom.copy_attributes(
                                                v[i].geom,
                                                [
                                                    "is_clean",
                                                    "latitude",
                                                    "longitude",
                                                    "source",
                                                ],
                                            )
                                            up = True
                                        elif (
                                            not v[i].geom.has_geom()
                                            and m.geom.has_geom()
                                        ):
                                            v[i].geom.shape = m.geom.shape
                                            v[i].geom.copy_attributes(
                                                m.geom,
                                                [
                                                    "is_clean",
                                                    "latitude",
                                                    "longitude",
                                                    "source",
                                                ],
                                            )
                                            up = True
                                        break
                                    else:
                                        nex = m.next_coms
                                elif m.mod == 21 and m.prev_coms != {k}:
                                    break
                                elif m.mod in [31, 32, 33]:
                                    lc = [
                                        *lc,
                                        *[n for n in m.prev_coms if n not in lc],
                                    ]
                                    nex = m.next_coms
                                elif m.mod in [20, 30, 41, 50, 70]:
                                    break
        return up

    def merge_geom(self, geoms: dict) -> bool:
        """Fusion de geometries"""
        up = False
        for k, v in geoms.items():
            # maj inutile
            if k.has_geom():
                continue

            # fusion
            if all([i.has_geom() for i in v]):
                up = True
                new_shape = unary_union([g.shape for g in v])
                k.shape = new_shape
                k.clean_geom()
                k.set_attribute("source", "Histocom")

        return up

    def geom_by_difference(self, geoms: dict) -> bool:
        """Geometrie par soustraction d'autres"""
        up = False
        for k, v in geoms.items():
            c = Counter([i.has_geom() for i in v])
            if k.has_geom() and c[False] == 1:
                geos = [g for g in v if g.has_geom()]
                new = [g for g in v if not g.has_geom()][0]
                new.shape = k.shape
                up = True
                for g in geos:
                    clip = g
                    new.shape = new.shape.difference(clip.shape)

                new.clean_geom()
                new.set_attribute("source", "Histocom")

        return up

    def geom_after_insee(self, paths: list):
        """Liste des communes apres INSEE"""
        for i, p in enumerate(paths):
            print("Après historique INSEE", p)
            dirdate = self.get_date(p)
            if dirdate <= self.insee_period[1]:
                paths_after = paths[: i + 1]
                paths_after.reverse()
                break

        current_coms = {}
        # verifier si "id" est renseigne pour les communes actuelles et creer une liste
        for k, v in self.communes.items():
            last_com = v[-1]
            if last_com.period_end is None:
                com_id = last_com.geom.attribute("id")
                if com_id is None:
                    raise ValueError(f"Missing id for com: {k}")
                else:
                    current_coms[com_id] = k

        # ajouter les communes qui n'ont jamais change
        last_insee_geom_path = paths_after.pop(0)
        geodf = gam.read_file(join(last_insee_geom_path, "COMMUNE.shp"))
        source = self.read_source(last_insee_geom_path)
        chfdf = None
        if isfile(join(last_insee_geom_path, "CHEF_LIEU.shp")):
            chfdf = gam.read_file(join(last_insee_geom_path, "CHEF_LIEU.shp"))

        if "ID" not in geodf:
            raise ValueError(f"There is no 'ID' field in: {last_insee_geom_path}")

        geodf = geodf.loc[lambda df: df["STATUT"] != "Arrondissement municipal"]
        for t in geodf.itertuples():
            if t.INSEE_COM not in self.communes:
                new_com = ComObj(
                    admin_code=t.INSEE_COM,
                    mod=99,
                    name=t.NOM_COM,
                    period_start=self.date_start,
                )
                new_geo = new_com.geom
                new_geo.set_attribute("id", t.ID)
                current_coms[t.ID] = t.INSEE_COM
                new_geo.set_attribute("source", source)
                if "POPULATION" in dir(t):
                    new_com.population = t.POPULATION

                new_geo.shape = t.geometry
                x_wgs, y_wgs = None, None
                if "X_CHF_LIEU" in dir(t):
                    project = Transformer.from_crs("EPSG:2154", "EPSG:4326")
                    x_l93, y_l93 = t.X_CHF_LIEU, t.Y_CHF_LIEU
                    if y_l93 < 1000000:
                        factor = 6 - int(log(y_l93, 10))
                        x_l93 = x_l93 * 10 ** factor
                        y_l93 = y_l93 * 10 ** factor

                    y_wgs, x_wgs = project.transform(x_l93, y_l93)
                elif chfdf is not None:
                    c2 = chfdf.loc[lambda df: df["INSEE_COM"] == t.INSEE_COM]
                    if c2["INSEE_COM"].count() == 1:
                        x_wgs, y_wgs = c2.iloc[0].geometry.x, c2.iloc[0].geometry.y

                new_geo.set_attribute("latitude", y_wgs)
                new_geo.set_attribute("longitude", x_wgs)

                self.communes[t.INSEE_COM] = [new_com]

        # faire les variations
        for gpath in paths_after:
            geodf = gam.read_file(join(gpath, "COMMUNE.shp"))
            source = self.read_source(gpath)
            chfdf = None
            if isfile(join(gpath, "CHEF_LIEU.shp")):
                chfdf = gam.read_file(join(gpath, "CHEF_LIEU.shp"))

            if "ID" not in geodf:
                raise ValueError(f"There is no 'ID' field in: {gpath}")

            geodf = geodf.loc[lambda df: df["STATUT"] != "Arrondissement municipal"]
            geodate = self.get_date(gpath)
            new_ids = list(geodf["ID"])
            new_coms = [i for i in new_ids if i not in current_coms]
            ending_coms = [i for i in current_coms.keys() if i not in new_ids]
            for c_id in ending_coms:
                c_insee = current_coms[c_id]
                self.communes[c_insee][-1].period_end = geodate - timedelta(days=1)

            for c_id in new_coms:
                record = geodf.loc[lambda df: df["ID"] == c_id].iloc[0]
                new_com = ComObj(
                    admin_code=record.INSEE_COM,
                    mod=99,
                    name=record.NOM_COM,
                    period_start=geodate,
                )
                new_geo = new_com.geom
                new_geo.set_attribute("id", record.ID)
                new_geo.set_attribute("source", source)
                new_geo.shape = record.geometry
                x_wgs, y_wgs = None, None
                if "X_CHF_LIEU" in record:
                    project = Transformer.from_crs("EPSG:2154", "EPSG:4326")
                    x_l93, y_l93 = record.X_CHF_LIEU, record.Y_CHF_LIEU
                    if y_l93 < 1000000:
                        factor = 6 - int(log(y_l93, 10))
                        x_l93 = x_l93 * 10 ** factor
                        y_l93 = y_l93 * 10 ** factor

                    y_wgs, x_wgs = project.transform(x_l93, y_l93)
                elif chfdf is not None:
                    c2 = chfdf.loc[lambda df: df["INSEE_COM"] == record.INSEE_COM]
                    if c2["INSEE_COM"].count() == 1:
                        x_wgs, y_wgs = c2.iloc[0].geometry.x, c2.iloc[0].geometry.y

                if "POPULATION" in record:
                    new_com.population = float(record.POPULATION)

                new_geo.set_attribute("latitude", y_wgs)
                new_geo.set_attribute("longitude", x_wgs)

                if record.INSEE_COM in self.communes:
                    self.communes[record.INSEE_COM].append(new_com)
                else:
                    self.communes[record.INSEE_COM] = [new_com]

    def update_coords(self, geo_to_update, insee_codes, period, type_period):
        for code in insee_codes:
            com_mod = self.get_mod(code, period, type_period)
            if com_mod is None:
                print(insee_codes, code, period, type_period)
                raise ValueError("è_é !")
            if com_mod.geom.attribute("longitude") is not None:
                geo_to_update.copy_attributes(com_mod.geom, ["longitude", "latitude"])
                return True

        return False

    def update_latlon(self):
        up = True
        while up:
            up = False
            for k, v in self.communes.items():
                for i, m in enumerate(v):
                    if m.geom.attribute("longitude") is None:
                        # a une geometrie
                        if m.geom.has_geom():
                            repr_pt = m.geom.shape.representative_point()
                            m.geom.set_attribute("longitude", repr_pt.x)
                            m.geom.set_attribute("latitude", repr_pt.y)
                            up = True
                        # derniere avec meme code insee
                        elif i + 1 < len(v):
                            if v[-1].geom.attribute("longitude") is not None:
                                m.geom.copy_attributes(
                                    v[-1].geom, ["longitude", "latitude"]
                                )
                                up = True
                        # recherche dans les suivantes
                        elif len(m.next_coms) > 0:
                            m1 = self.update_coords(
                                m.geom,
                                m.next_coms,
                                m.period_end + timedelta(days=1),
                                ComObj.P_START,
                            )
                            up = True if m1 else up
                        # recherche dans les precedentes
                        elif len(m.prev_coms) > 0:
                            m1 = self.update_coords(
                                m.geom,
                                m.prev_coms,
                                m.period_start - timedelta(days=1),
                                ComObj.P_END,
                            )
                            up = True if m1 else up

    def df(self):
        data = {
            "insee_code": [],
            "mod": [],
            "name": [],
            "period_start": [],
            "period_end": [],
            "previous": [],
            "next": [],
            "population": [],
            "geom": [],
        }

        for k, v in self.communes.items():
            for m in v:
                data["insee_code"].append(k)
                data["mod"].append(m.mod)
                data["name"].append(m.name)
                data["period_start"].append(m.period_start)
                data["period_end"].append(m.period_end)
                data["previous"].append(m.prev_coms)
                data["next"].append(m.next_coms)
                data["population"].append(m.population)
                data["geom"].append(m.geom)

        df = am.DataFrame(data)
        return df

    def tz_to_utc(self, insee_code, date_time):
        if date_time is None:
            return None

        if insee_code[:3] in ["971", "972", "973", "974", "976"]:
            tz = self.timezones[insee_code[:3]]
        else:
            tz = self.timezones["paris"]

        tz_f = "%Y-%m-%dT%H:%M:%S.%f"
        dt_tz = tz.localize(date_time)
        dt_utc = dt_tz.astimezone(self.timezones["utc"])
        return dt_utc.strftime(tz_f)

    def export_csv(self):
        with open(join(os.getcwd(), "export_notebook.csv"), "w") as f:
            f.write(
                '"administrative_code","code_ref","name","period_start","period_end"\n'
            )
            for k, v in self.communes.items():
                for m in v:
                    period_start = datetime.combine(m.period_start, datetime.min.time())
                    if m.period_end is None:
                        period_end = None
                    else:
                        period_end = datetime.combine(
                            m.period_end, datetime.min.time()
                        ) - timedelta(microseconds=1)

                    f.write(
                        f'"{k}",'
                        f'"{",".join(m.next_coms)}",'
                        f'"{m.name}",'
                        f'"{self.tz_to_utc(k, period_start)}",'
                        f'"{self.tz_to_utc(k, period_end)}"\n'
                    )

    def export_geojson(self, path: str):
        data = {
            "country_iso2": [],
            "source": [],
            "source_version": [],
            "source_id": [],
            "origin_source": [],
            "origin_source_id": [],
            "administrative_code_type": [],
            "administrative_code": [],
            "population": [],
            "district": [],
            "name": [],
            "name_fr": [],
            "period_start": [],
            "period_end": [],
            "longitude": [],
            "latitude": [],
            "geometry": [],
        }

        source_version = str(datetime.now().date())

        for k, v in self.communes.items():
            for m in v:
                print(m.name, m.period_start)
                period_start = datetime.combine(m.period_start, datetime.min.time())
                if m.period_end is None:
                    period_end = None
                else:
                    period_end = datetime.combine(
                        m.period_end, datetime.min.time()
                    ) - timedelta(microseconds=1)

                data["country_iso2"].append("FR")
                data["source"].append("Histocom")
                data["source_version"].append(source_version)
                data["source_id"].append(f"{k}_{m.period_start.strftime('%Y_%m_%d')}")
                data["origin_source"].append(m.geom.attribute("source"))
                data["origin_source_id"].append(str(m.geom.attribute("id")))
                data["administrative_code_type"].append("INSEE")
                data["administrative_code"].append(k)
                data["population"].append(m.population)
                data["district"].append(False)
                data["name"].append(m.name)
                data["name_fr"].append(m.name)
                data["period_start"].append(self.tz_to_utc(k, period_start))
                data["period_end"].append(self.tz_to_utc(k, period_end))
                data["longitude"].append(m.geom.attribute("longitude"))
                data["latitude"].append(m.geom.attribute("latitude"))
                if m.geom.shape is None:
                    data["geometry"].append(
                        Point(
                            m.geom.attribute("longitude"),
                            m.geom.attribute("latitude"),
                        )
                    )
                elif m.geom.shape.geom_type == "Polygon":
                    data["geometry"].append(MultiPolygon([m.geom.shape]))
                else:
                    data["geometry"].append(m.geom.shape)

        gdf = gam.GeoDataFrame(data)
        # outre-mer
        gdf_outre_mer = gam.read_file(f"{path}/collectivites_outre_mer.geojson")
        # assemblage
        gdf = gdf.append(gdf_outre_mer)
        # export
        gdf.to_file(
            f"{path}/histocom_{source_version.replace('-', '')}.geojson", driver="GeoJSON"
        )

    def count(self):
        i = 0
        for v in self.communes.values():
            i += len(v)
        return i


workdir = "coms/data"
dir_insee = join(workdir, "INSEE")
dir_geom = join(workdir, "GEOM")

if isdir(dir_insee):
    insee_files = [f for f in os.scandir(dir_insee) if isfile(f)]
    if len(insee_files) > 0:
        insee = insee_files[0].path
    else:
        raise ValueError("Wrong path")

histo = Histocom()
print("Build INSEE data")
histo.build_insee(insee)
print("Build geometry data")
histo.build_geom(dir_geom)
print("Export to geoJSON")
histo.export_geojson("coms")
