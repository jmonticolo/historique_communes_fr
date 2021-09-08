#!/bin/bash

WDIR=data
GET_DATA=1
F_CHL="CHEF_LIEU"
F_COM="COMMUNE"

AE_L93_UTM="\
ADMIN-EXPRESS_2-0__SHP__FRA_2019-03-14 \
ADMIN-EXPRESS_2-0__SHP__FRA_2019-02-15 \
ADMIN-EXPRESS_2-0__SHP__FRA_2019-01-18\
"

AE_WGS="\
ADMIN-EXPRESS_2-3__SHP__FRA_WM_2020-05-18 \
ADMIN-EXPRESS_2-3__SHP__FRA_WM_2020-03-16 \
ADMIN-EXPRESS_2-2__SHP__FRA_WM_2020-02-24 \
ADMIN-EXPRESS_2-2__SHP__FRA_WM_2020-01-16 \
ADMIN-EXPRESS_2-1__SHP__FRA_WM_2019-12-16 \
ADMIN-EXPRESS_2-1__SHP__FRA_WM_2019-11-15 \
ADMIN-EXPRESS_2-1__SHP__FRA_WM_2019-10-15 \
ADMIN-EXPRESS_2-1__SHP__FRA_WM_2019-09-16 \
ADMIN-EXPRESS_2-1__SHP__FRA_WM_2019-08-19 \
ADMIN-EXPRESS_2-0__SHP__FRA_WM_2019-07-16 \
ADMIN-EXPRESS_2-0__SHP__FRA_WM_2019-06-18 \
ADMIN-EXPRESS_2-0__SHP__FRA_WM_2019-05-15\
"

AE_2020="ADMIN-EXPRESS-COG_2-1__SHP__FRA_WGS84_2020-03-25"
AE_2019="ADMIN-EXPRESS-COG_2-0__SHP__FRA_WGS84G_2019-09-24"
AE_2018="2018-05-04\$ADMIN-EXPRESS-COG_1-1__SHP__FRA_2018-04-03/file/ADMIN-EXPRESS-COG_1-1__SHP__FRA_2018-04-03"
AE_2017="2017-07-07\$ADMIN-EXPRESS-COG_1-0__SHP__FRA_2017-06-19/file/ADMIN-EXPRESS-COG_1-0__SHP__FRA_2017-06-19"

GEOFLA_2016="GEOFLA_THEME-COMMUNE_2016\$GEOFLA_2-2_*_2016-06-28/file/GEOFLA_2-2_*_2016-06-28"
GEOFLA_2015="GEOFLA_THEME-COMMUNE_2015_2\$GEOFLA_2-1_*_2015-12-01/file/GEOFLA_2-1_*_2015-12-01"
GEOFLA_REGIONS="\
COMMUNE_SHP_LAMB93_FXX \
COMMUNE_SHP_UTM20W84GUAD_D971 \
COMMUNE_SHP_UTM20W84MART_D972 \
COMMUNE_SHP_UTM22RGFG95_D973 \
COMMUNE_SHP_RGR92UTM40S_D974 \
COMMUNE_SHP_RGM04UTM38S_D976\
"

GEOFLA_2014="GEOFLA_THEME-COMMUNE_2014_GEOFLA_2-0_COMMUNE_SHP_LAMB93_FXX_2015-03-23/file/GEOFLA_2-0_COMMUNE_SHP_LAMB93_FXX_2015-03-23"
GEOFLA_2013="GEOFLA_THEME-COMMUNE_2013_GEOFLA_1-1_SHP_LAMB93_FR-ED131/file/GEOFLA_1-1_SHP_LAMB93_FR-ED131"
GEOFLA_2012="GEOFLA_THEME-COMMUNE_2012_GEOFLA_1-1_SHP_LAMB93_000_2013-01-15/file/GEOFLA_1-1_SHP_LAMB93_000_2013-01-15"
GEOFLA_2011="GEOFLA_THEME-COMMUNE_2011_GEOFLA_1-1_SHP_LAMB93_FR-ED111/file/GEOFLA_1-1_SHP_LAMB93_FR-ED111"

mkdir $WDIR
mkdir $WDIR/dl
mkdir $WDIR/GEOM
mkdir $WDIR/INSEE

if [ $GET_DATA == 1 ]
  then
    # GET DATA
    ## INSEE
    URL="https://www.insee.fr/fr/statistiques/fichier/4316069/mvtcommune2020-csv.zip"
    wget -P ./$WDIR/dl/ $URL
    7z e ./$WDIR/dl/mvtcommune2020-csv.zip -o$WDIR/INSEE/

    ## GEOFLA 2011-2014
    URL="https://wxs.ign.fr/oikr5jryiph0iwhw36053ptm/telechargement/inspire"
    for file in $GEOFLA_2011 $GEOFLA_2012 $GEOFLA_2013 $GEOFLA_2014
    do
        wget -P ./$WDIR/dl/ "$URL/$file.7z"
    done

    ## GEOFLA 2015-2016
    for version in $GEOFLA_2015 $GEOFLA_2016
    do
        YEAR=$(echo $version | rev | cut -d _ -f 1 | rev | cut -d - -f 1)
        for file in $GEOFLA_REGIONS
        do
            zone=$(echo $file | rev | cut -d _ -f 1 | rev)
            if [ $YEAR == "2016" ] && [ $zone == "D976" ]
                then
                    file="_SHP_RGM04UTM38S_D976"
                    wget -O "./$WDIR/dl/GEOFLA_2-2_COMMUNE"$file"_2016-06-28.7z" "$URL/${version//\*/$file}.7z"
                else
                    wget -P ./$WDIR/dl/ "$URL/${version//\*/$file}.7z"
            fi
        done
    done

    ## ADMIN-EXPRESS 2017-2018
    URL="https://wxs.ign.fr/x02uy2aiwjo9bm8ce5plwqmr/telechargement/prepackage"
    for file in $AE_2017 $AE_2018
    do
        wget -P ./$WDIR/dl/ "$URL/ADMINEXPRESS-COG-PACK_$file.7z"
    done

    ## ADMIN-EXPRESS >= 2019
    URL="ftp://Admin_Express_ext:Dahnoh0eigheeFok@ftp3.ign.fr"
    for file in $AE_2019 $AE_2020 $AE_L93_UTM $AE_WGS
    do
        wget -P ./$WDIR/dl/ "$URL/$file.7z.001"
    done
fi

# EXTRACT DATA
## IGN L93_UTM MONTH
for file in $AE_L93_UTM
do
    MONTH=$(echo $file | rev | cut -d _ -f 1 | rev | cut -d - -f 2)
    YEAR=$(echo $file | rev | cut -d _ -f 1 | rev | cut -d - -f 1)
    DIR="./$WDIR/GEOM/$YEAR-$MONTH"
    SOURCE=$(echo $file | cut -d _ -f 1)
    echo $SOURCE > $DIR/source.txt
    for shp in $F_CHL $F_COM
    do
        for zone in {FR,D971,D972,D973,D974,D976}
            do
                7z e ./$WDIR/dl/$file.7z.001 -o$DIR/$zone/ ADMIN*/*/1_*/ADE*_$zone/$shp.*
                chmod 755 $DIR/$zone/*
                if [ -f "$DIR/$shp.shp" ]
                    then
                        ogr2ogr -f "ESRI Shapefile" $DIR/$shp.shp $DIR/$zone/$shp.shp -update -append -t_srs EPSG:4326 -progress
                    else
                        ogr2ogr -f "ESRI Shapefile" $DIR/$shp.shp $DIR/$zone/$shp.shp -t_srs EPSG:4326 -progress
                fi
            done
    done
done

## IGN WGS84 MONTH
for file in $AE_WGS
do
    MONTH=$(echo $file | rev | cut -d _ -f 1 | rev | cut -d - -f 2)
    YEAR=$(echo $file | rev | cut -d _ -f 1 | rev | cut -d - -f 1)
    DIR="./$WDIR/GEOM/$YEAR-$MONTH/"
    SOURCE=$(echo $file | cut -d _ -f 1)
    echo $SOURCE > $DIR/source.txt
    for shp in $F_CHL $F_COM
    do
        7z e ./$WDIR/dl/$file.7z.001 -o$DIR ADMIN*/*/1_*/ADE*/$shp.*
        chmod 755 $DIR*
    done
done

## IGN WGS84 YEAR
for file in $AE_2019 $AE_2020
do
    YEAR=$(echo $file | rev | cut -d _ -f 1 | rev | cut -d - -f 1)
    DIR="./$WDIR/GEOM/$YEAR/"
    SOURCE=$(echo $file | cut -d _ -f 1)
    echo $SOURCE > $DIR/source.txt
    for shp in $F_CHL $F_COM
    do
        7z e ./$WDIR/dl/$file.7z.001 -o$DIR ADMIN*/*/1_*/ADE*/$shp.*
        chmod 755 $DIR*
    done
done

## IGN L93_UTM YEAR
for rep in $AE_2017 $AE_2018
do
    file=$(echo $rep | rev | cut -d / -f 1 | rev)
    YEAR=$(echo $file | rev | cut -d _ -f 1 | rev | cut -d - -f 1)
    DIR="./$WDIR/GEOM/$YEAR"
    SOURCE=$(echo $file | cut -d _ -f 1)
    echo $SOURCE > $DIR/source.txt
    for shp in $F_CHL $F_COM
    do
        for zone in {FR,D971,D972,D973,D974,D976}
            do
                7z e ./$WDIR/dl/$file.7z -o$DIR/$zone/ ADMIN*/*/1_*/ADE*_$zone/$shp.*
                chmod 755 $DIR/$zone/*
                if [ -f "$DIR/$shp.shp" ]
                    then
                        ogr2ogr -f "ESRI Shapefile" $DIR/$shp.shp $DIR/$zone/$shp.shp -update -append -t_srs EPSG:4326 -progress
                    else
                        ogr2ogr -f "ESRI Shapefile" $DIR/$shp.shp $DIR/$zone/$shp.shp -t_srs EPSG:4326 -progress
                fi
            done
    done
done

## GEOFLA 2011-2014
for version in "GEOFLA_2011" "GEOFLA_2012" "GEOFLA_2013" "GEOFLA_2014"
do
    rep=${!version}
    file=$(echo $rep | rev | cut -d / -f 1 | rev)
    YEAR=$(echo $version | rev | cut -d _ -f 1 | rev)
    DIR="./$WDIR/GEOM/$YEAR"
    SOURCE=$(echo $file | cut -d _ -f 1)
    echo $SOURCE > $DIR/source.txt
    for shp in $F_COM
    do
        7z e ./$WDIR/dl/$file.7z -o$DIR/FXX/ GEOFLA*/*/1_*/GEOFLA*/*/$shp.*
        chmod 755 $DIR/FXX/*
        ogr2ogr -f "ESRI Shapefile" $DIR/$shp.shp $DIR/FXX/$shp.SHP -t_srs EPSG:4326 -progress
    done
done

## GEOFLA 2015-2016
for version in $GEOFLA_2015 $GEOFLA_2016
do
    file_template=$(echo $version | rev | cut -d / -f 1 | rev)
    for reg in $GEOFLA_REGIONS
    do
        file=${file_template//\*/$reg}
        zone=$(echo $reg | rev | cut -d _ -f 1 | rev)
        YEAR=$(echo $file | rev | cut -d _ -f 1 | rev | cut -d - -f 1)
        DIR="./$WDIR/GEOM/$YEAR"
        SOURCE=$(echo $file | cut -d _ -f 1)
        echo $SOURCE > $DIR/source.txt
        for shp in $F_COM
        do
            7z e ./$WDIR/dl/$file.7z -o$DIR/$zone/ GEOFLA*/*/1_*/GEOFLA*/*/$shp.*
            chmod 755 $DIR/$zone/*
            if [ -f "$DIR/$shp.shp" ]
                then
                    ogr2ogr -f "ESRI Shapefile" $DIR/$shp.shp $DIR/$zone/$shp.shp -update -append -t_srs EPSG:4326 -progress
                else
                    ogr2ogr -f "ESRI Shapefile" $DIR/$shp.shp $DIR/$zone/$shp.shp -t_srs EPSG:4326 -progress
            fi
        done
    done
done
