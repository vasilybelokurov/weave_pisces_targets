
query = 'SELECT gs.source_id, xm.ra,xm.dec,xm.j_m,xm.h_m,xm.k_m,gs.ebv,xm.pmra,xm.pmdec,xm.pmra_error,xm.pmdec_error,xm.phot_g_mean_mag,xm.bp_rp,xm.parallax,xm.parallax_error,xm.phot_bp_rp_excess_factor,xm.phot_g_n_obs FROM gaia_edr3_aux.gaia_source_2mass_xm as xm, gaia_edr3.gaia_source as gs where xm.source_id=gs.source_id and xm.parallax <0.5 AND abs(xm.b) > 5 and xm.pmra^2+xm.pmdec^2 < 25;'

d = pgsql_query(query)

save, filename = '~/data/gaia/edr3/low_parallax05_b5_pm5_2mass_ebv.idl', d


end
