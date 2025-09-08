# creating the gaia HDF5 


* Footprint notes
We assume >-10 and |b|>10


For observing conditions for both pointed and shared surveys we had OBSTEMP=DACEB : seeing<1", sky transparency>0.8, airmass<1.5, no constraint on moon distance (this should already be handled by sky brightness constraint), sky brightness > 21.5 mag/sq arcsec in V.

The nominal PROGTEMP for the shared is 11331+ (MOS LR, 3x20min, 1x binning, chained OB).
However the current plan for the binning is that we submit both binned and unbinned OBs, and catalogues, and then retract the ones that we don't want once the decision is made. So we also need each row duplicated with at least 11332+ for the binned option. This is obviously going to impact the catalogue size...

For the pointed we always want 11331.X where X is the number of times that we want that field observed in the semester. We need to discuss the plan for X, since it won't just be the number of times we want to observe that target because then we get no timebase control for binary detection.

We need to discuss the OBSTEMP and PROGTEMP values with WL and WQ since we need make sure we choose the same values otherwise I think I can't merge the catalogues together to make OBs.


***********



https://gea.esac.esa.int/archive/documentation/GEDR3/Catalogue_consolidation/chap_crossmatch/sec_crossmatch_externalCat/ssec_panstarrs.html






##### New selections

This is the selection of good quality panstarrs objects (see gaia edr3 page on xmatches)

create table  panstarrs_sub3 as select o.objid, ramean,decmean, 
        gmeanpsfmag, rmeanpsfmag, imeanpsfmag, zmeanpsfmag, ymeanpsfmag, gmeankronmag,
        rmeankronmag,imeankronmag,zmeankronmag,ymeankronmag 
        from panstarrs_dr1.objectthin as o, panstarrs_dr1.meanobject as mo  where ndetections!=1 and  
        (objInfoFlag&33554432>0) and  (objInfoFlag&524288)=0 and (objInfoFlag&1048576)=0 and o.objid=mo.objid;

# This is crossmatchign 
create table gaia_dr3_ps1_xm_2 as
        select objid,
        gmeanpsfmag, rmeanpsfmag, imeanpsfmag, zmeanpsfmag, ymeanpsfmag, 
        gmeankronmag,rmeankronmag,imeankronmag,zmeankronmag,ymeankronmag,
       ramean as p_ra ,decmean as p_dec, g.*
       from panstarrs_sub3 as ps left join
            gaia_edr3_aux.panstarrs1bestneighbour  as bs
        on (bs.original_ext_source_id=ps.objid)
        left join gaia_edr3.gaia_source as g on
	   (bs.source_id= g.source_id);


IFS='' read -r -d '' QUERY <<"EOF"
with gwise as (
select source_id,phot_g_mean_mag,ra,dec from gaiadr3_wise 
where (w1mpro-w2mpro-2*sqrt(w1sigmpro^2+w2sigmpro^2))> 
greatest(0.5,0.5-(phot_g_mean_mag-2*ebv-w1mpro-2.75)) 
) 
select  objid, coalesce(source_id,-1) as source_id, 
   coalesce(ra,p_ra) as ra , coalesce(dec,p_dec) as dec,
   phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,pmra, 
   pmdec,parallax, parallax_error, 
   pmra_error,pmdec_error,pmra_pmdec_corr,
   gmeanpsfmag,rmeanpsfmag, 
   gmeankronmag,rmeankronmag,astrometric_excess_noise, ruwe, 
   healpix_ang2ipix_nest(128,coalesce(ra,p_ra),coalesce(dec,p_dec)) as hpx128,
   exists (select 1 from gwise as g where gp.source_id=g.source_id) as qso_flag
   from 
    gaia_dr3_ps1_xm_2 as gp where coalesce(dec,p_dec)>-10 and 
   abs(fk52galb(coalesce(ra,p_ra),coalesce(dec,p_dec))) > 10  
   order by healpix_ang2ipix_nest(128,coalesce(ra,p_ra),coalesce(dec,p_dec))
EOF

python ~/astrolibpy/my_utils/pg2hdf5.py --host 'cappc127.ast.cam.ac.uk' "$QUERY" ../data_repo/gaiaps1_edr3_220425.h5

gaia_split ../data_repo/gaiaps1_edr3_220425.h5  ../data_repo/gaiaps1edr3_220425_hpx/ gaiaps1

IFS='' read -r -d '' QUERY1 <<"EOF"
with x as 
(
   select source_id,ra,dec,parallax,parallax_error,
   pmra,pmdec,pmra_error,pmdec_error, 0 as pmra_pmdec_corr,
   phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,
   sloan_g,sloan_r,popid, log_grav from gaiaedr3mock.main  as g,
   gaiaedr3mock.parsec_props as p where (b>10 or b<-10) and dec>-10 
   and g.index_parsec=p.parsec_index
) 
select row_number() over() as objid,source_id, 
ra, dec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
pmra+sqrt(-2*ln(random()))*cos(2*3.141592*random()) * pmra_error as pmra,
pmdec+sqrt(-2*ln(random()))*cos(2*3.141592*random()) * pmdec_error as pmdec,
parallax+sqrt(-2*ln(random()))*cos(2*3.141592*random()) * parallax_error as parallax,
pmra_pmdec_corr,
parallax_error, pmra_error,pmdec_error,
sloan_g-0.011-0.125*(sloan_g-sloan_r)-0.015*(sloan_g-sloan_r)^2 + 5*log(1e3/parallax)-5 as  gmeanpsfmag,
sloan_r+0.001-0.006*(sloan_g-sloan_r)-0.002*(sloan_g-sloan_r)^2 + 5*log(1e3/parallax)-5 as rmeanpsfmag ,
sloan_g-0.011-0.125*(sloan_g-sloan_r)-0.015*(sloan_g-sloan_r)^2 + 5*log(1e3/parallax)-5 as  gmeankronmag,
sloan_r+0.001-0.006*(sloan_g-sloan_r)-0.002*(sloan_g-sloan_r)^2 + 5*log(1e3/parallax)-5 as rmeankronmag,
0. as astrometric_excess_noise, 1. as ruwe,
healpix_ang2ipix_nest(128,ra,dec) as hpx128,
popid, log_grav, 5*log(1e3/parallax)-5 as dm from x
order by healpix_ang2ipix_nest(128,ra,dec)
EOF

python ~/astrolibpy/my_utils/pg2hdf5.py --host 'cappc127.ast.cam.ac.uk' "$QUERY1" ../data_repo/gaiaedr3_mock_220425.h5
gaia_split ../data_repo/gaiaedr3_mock_220425.h5 ../data_repo/gaiaedr3_mock/ gaiaps1



************ 
generating allsky/mock 


weave_make_blue_red --input_prefix /data2/weave/catalog_selection/data_repo/gaiaedr3_mock/ --output_prefix ${ROOT}/data_repo/internal_cats/testing/ --version 200000 --no_extinction_correct --fields_fits ${ROOT}/data_repo/footprints/fake_NGP_footprint.fits



wisexGaia

create table gaiadr3_wise as with x as (select *,ebv from gaia_dr3.gaia_source as g where pmra is null or (not (abs(pmra)> 3*pmra_error) and not (abs(pmdec)>3*pmdec_error) and not (parallax> 3*parallax_error))), y as materialized ( select source_id,ra,dec,phot_g_mean_mag,w1mpro,w2mpro,w1sigmpro,w2sigmpro,ebv from x, lateral (select w1mpro, w2mpro,w1sigmpro,w2sigmpro from catwise_202003.main as w where q3c_join(x.ra,x.dec,w.ra,w.dec,1./3600) order by q3c_dist(x.ra,x.dec,w.ra,w.dec) asc limit 1) as aa) select * from y where w1mpro is not null;

