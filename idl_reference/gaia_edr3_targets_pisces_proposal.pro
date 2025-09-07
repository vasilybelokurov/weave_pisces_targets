
; SETUP

tarx = [-45, 15]
tary = [-20, 0]

rpm_lim = [0, 8.]+10  ; [0, 6.]+10
par_cut = 0.25 ; 0.1
mag_lim = [13, 18.5] ; [14, 18]
col_lim = [1., 4] ;  [1.25, 4]
dist_cut = 50

rpm_lim2 = [0, 8]+10
par_cut2 = 0.2
mag_lim2 = [14, 20.]
col_lim2 = [1.1, 4]
dist_cut2 = 30


dmagr = [15, 21]
distr = [10, 110]
lr = [-180, 180]

rar = [-180, 180]
decr = [-80, 80]
nra = 70
ndec = 35

rlr = reverse(lr)
rrar = reverse(rar)

perc = [30, 90]

fwhm = 1.

nra_w = 80
ndec_w = 60
fwhm_w = 1.3
distcut_w = [65, 100]

; -----------------

; 2MASS
coeff = 1.71
off = 0.68
doffm = 0.2
doffp = 0.15
hklim = [0.04, 10]

; Distances
dist_lmc = 49.97
dm_lmc = 5*alog10(dist_lmc*1e+3)-5
dm_smc = 5*alog10(62.1*1e+3)-5
dm_sgr = 5*alog10(28*1e+3)-5

col_nodes_lmc = [0.97, 1.27, 1.5, 1.8, 2.5, 3.5]
mag_nodes_lmc = [20,  17.3, 16.3, 15.7, 15.3, 15.4]

;; col_nodes_smc = [0.97, 1.27, 1.5, 1.8, 2.5, 3.5]
;; mag_nodes_smc = [19.7,  17.15, 16.3, 15.7, 15.6, 15.8]

col_nodes_smc = [0.97, 1.27,  1.5,   1.8,   2.5, 3.5]
mag_nodes_smc = [19.7, 17.25, 16.4, 15.8, 15.62, 15.8]

col_nodes_sgr = [0.97, 1.27, 1.5,    1.8,  2.1, 2.5, 3.5]
mag_nodes_sgr = [19.5, 16.48, 15.37, 14.8, 14.6, 14.5, 14.5]

;GOTO, plot
;GOTO, loaded

restore, '~/data/gaia/edr3/low_parallax05_b5_pm5_2mass_ebv.idl'

; Cut known satellites
in = cutradec(d.ra, d.dec, rad = 0.5, out = out, /comp)


; Extinction
A0 = 3.1*d.ebv
kg = 0.9761-0.1704*d.bp_rp
kbp = 1.1517-0.0871*d.bp_rp
krp = 0.6104-0.0170*d.bp_rp
AG = kg*A0
ABP = kbp*A0
ARP = krp*A0
mg = d.PHOT_G_MEAN_MAG-ag
col = d.bp_rp-abp+arp
kg = 0 & kbp = 0 & krp = 0 & ag = 0 & abp = 0 & arp = 0 & a0 = 0

jmag = d.j_m-0.72*d.ebv
hmag = d.h_m-0.46*d.ebv
kmag = d.k_m-0.306*d.ebv

; Distances
dmag_lmc = interpol(mag_nodes_lmc, col_nodes_lmc, col, /spline)-dm_lmc
dmag_smc = interpol(mag_nodes_smc, col_nodes_smc, col, /spline)-dm_smc
dmag_sgr = interpol(mag_nodes_sgr, col_nodes_sgr, col, /spline)-dm_sgr

dmag_lmc = mg-dmag_lmc
dmag_smc = mg-dmag_smc
dmag_sgr = mg-dmag_sgr

dmag = 0.5*(dmag_sgr+dmag_smc)

; Galactic
euler, d.ra, d.dec, l, b, 1
l = w360to180(l)

; Sgr
ori_sgr = [283.750, -30.483]        ; at Sgr's location
pole_sgr = [303.63, 59.58]
gc = sphcrd_trans(d.ra, d.dec, pole_sgr[0], pole_sgr[1], c10 = ori_sgr[0], c20 = ori_sgr[1])
xsgr = reform(gc[*, 0])
ysgr = reform(gc[*, 1])
gc = 0

; PM
pm = sqrt(d.pmra^2+d.pmdec^2)
;pm2 = sqrt(pml2^2+pmb2^2)

; RPM
rpm = mg+5.*alog10(pm)-1.47*sin(abs(b*!dtor))

; SGR simulation
dir_nolmc = '~/data/sims/eugene/Sgr+MW_noLMC/'
dir_lmc = '~/data/sims/eugene/Sgr+LMC+MW/'

restore, '~/data/sims/eugene/template.idl'
dlmc = read_ascii(dir_lmc+'stars.txt', template = template)
euler, dlmc.ra, dlmc.dec, llmc, blmc, 1
dnolmc = read_ascii(dir_nolmc+'stars.txt', template = template)
euler, dnolmc.ra, dnolmc.dec, lnolmc, bnolmc, 1

lalmc = dlmc.lambda
w = where(lalmc LT -180)
lalmc[w] = 360+lalmc[w]
w = where(lalmc GT 180)
lalmc[w] = lalmc[w]-360

lanolmc = dnolmc.lambda
w = where(lanolmc LT -180)
lanolmc[w] = 360+lanolmc[w]
w = where(lanolmc GT 180)
lanolmc[w] = lanolmc[w]-360

wlmc = where(dlmc.lambda LT -360)
wnolmc = where(dnolmc.lambda LT -360)

; Zaritsky
path = '~/Documents/Work/lists/zaritsky_h3.dat'
n = 15
dz = fltarr(n, 3)
openr, 1, path
readf, 1, dz
close, 1
ra_z = reform(dz[*, 0])
dec_z = reform(dz[*, 1])
v_z = reform(dz[*, 2])

; WAKE sim
restore, '~/data/sims/denis/final_stellar_halo_1.5e11Msun_LMC.idl'
euler, dw.l, dw.b, ra_w, dec_w, 2

restore, '~/data/sims/denis/lmc_pos.template'
lmco = read_ascii('~/data/sims/denis/lmc_pos.txt', template = template)
euler, lmco.l, lmco.b, ra_lmc, dec_lmc, 2


loaded:

; SEELCTION
wpar = where(d.parallax LT par_cut)
wmag = where(mg GT mag_lim[0] AND mg LT mag_lim[1])
wcol = where(col GT col_lim[0] AND col LT col_lim[1])
wrpm = where(rpm GT rpm_lim[0] AND rpm LT rpm_lim[1])

; targets
wpar2 = where(d.parallax LT par_cut2)
wmag2 = where(mg GT mag_lim2[0] AND mg LT mag_lim2[1])
wcol2 = where(col GT col_lim2[0] AND col LT col_lim2[1])
wrpm2 = where(rpm GT rpm_lim2[0] AND rpm LT rpm_lim2[1])

wxy = where(w360to180(d.ra) GT tarx[0] AND w360to180(d.ra) LT tarx[1] AND d.dec GT tary[0] AND d.dec LT tary[1], ntar_xy)

; 2MASS
wjhk = where(jmag-kmag GT coeff*(hmag-kmag)+off-doffm AND jmag-kmag LT coeff*(hmag-kmag)+off+doffp)
whk = where(hmag-kmag GT hklim[0] AND hmag-kmag LT hklim[1])
; combine
w2mass = isection([wjhk, whk], 2)

; Giant selection
wg = isection([wpar, wmag, wcol, wrpm], 4)


wg_tar = isection([wpar2, wmag2, wcol2, wrpm2], 4)
wg2 = isection([wpar, wmag, wcol, wjhk, wrpm], 5)
wsgr = where(abs(ysgr) LT 15)
wdist = where(10^(0.2*(dmag+5))*1e-3 GT dist_cut)

wg_sgr = isection([wg, wsgr], 2)
wg2_sgr = isection([wg2, wsgr], 2)

wtar = isection([wxy, wg_tar], 2, n = ntar)
wtar2 = isection([wxy, out, wg], 3, n = ntar2)
wtar3 = isection([wxy, out, wg, wdist], 4, n = ntar3)

ww = isection([out, wg, wdist], 3)

den_radec = hist2d(w360to180((d.ra)[ww]), (d.dec)[ww], rar, decr, nra, ndec, /xflip)

;; ww = isection([out, wg2, wdist], 3)
;; den_radec2 = hist2d(w360to180((d.ra)[ww]), (d.dec)[ww], rar, decr, nra, ndec, /xflip)

plot:

;; w = where(dw.dist GT distcut_w[0] AND dw.dist LT distcut_w[1] AND w360to180(ra_w) GT -70)
;; denw_radec = hist2d(w360to180(ra_w[w]), dec_w[w], rar, decr, nra_w, ndec_w, /xflip)


set_plot, 'ps'
DEVICE, filename = '~/work/plots/gaia_edr3_lmc/pisces_giants.ps', xsize = 9, ysize = 3, /inches, bits_per_pixel = 32, /color, xoffset = 0.5, /encapsulated

!p.font = 0
!p.multi = [0, 2, 1]
!p.charsize = 1.
!x.margin = [6.5, 0.5]
!y.margin = [3.5, 1]

plot, [0], [0], psym = 3, xr = rrar, xs = 1, yr = decr, ys = 1, xtickformat = '(a1)', ytickformat = '(a1)'
im = float(den_radec)
w0 = where(im GT 0)
im = filter_image(im, fwhm = fwhm, /all)
im = alog10(im)
mm = prank(im[w0], perc)
tvimage, 255-bytscl(im, min = mm[0], max = mm[1]), /overplot, /noint
axis_oplot, xr = rrar, yr = decr, xtitle = 'RA', ytitle = 'Dec'

cgLoadCT, 0
LoadCT, 31, FILE='~/idl/extra/fsc_brewer.tbl', bottom = 1, ncolors = 254

;; imw = float(denw_radec)
;; imw = filter_image(imw, fwhm = 1.3, /all)
;; w0 = where(imw GT 0)
;; mm = prank(imw[w0], perc)
;; xy = pix2d(rrar, decr, nra_w, ndec_w)
;; c_lev = prank(imw, [98.2, 99.2])
;; contour, imw, xy.xx, xy.yy, levels = c_lev, /overplot, c_color = 254

;; ww = where(dec_lmc LT -20)
;; oplot, w360to180(ra_lmc[ww]), dec_lmc[ww], color = 254

loadct, 0

; -------------------

plot, [0], [0], psym = 3, xr = rrar, xs = 1, yr = decr, ys = 1, xtickformat = '(a1)', ytickformat = '(a1)'
im = float(den_radec)
w0 = where(im GT 0)
im = filter_image(im, fwhm = fwhm, /all)
im = alog10(im)
mm = prank(im[w0], perc)
tvimage, 255-bytscl(im, min = mm[0], max = mm[1]), /overplot, /noint
axis_oplot, xr = rrar, yr = decr, xtitle = 'RA', ytitle = 'Dec'

oplot, rrar, [0, 0], color = 0, linestyle = 0

cgLoadCT, 0
LoadCT, 31, FILE='~/idl/extra/fsc_brewer.tbl', bottom = 1, ncolors = 254

; Sgr
oplot, w360to180((dnolmc.ra)[wnolmc]), (dnolmc.dec)[wnolmc], psym = sym(1),  symsize = 0.2, color = 154

; Zaritsky
oplot, w360to180(ra_z), dec_z, psym = sym(1), symsize = 0.45, color = 200

; Sgr
oplot, w360to180((dlmc.ra)[wlmc]), (dlmc.dec)[wlmc], psym = sym(1), symsize = 0.3, color = 60

; Target
oplot, tarx[[0, 1, 1, 0, 0]], tary[[0, 0, 1, 1, 0]], color = 254, thick = 2

xyouts, 49, 1, 'Pisces', color = 254, charsize = 0.7
xyouts, 75, -7, 'Overdensity', color = 254, charsize = 0.7

xyouts, 135, 22, 'Sgr Tr', charsize = 0.7, color = 255
xyouts, -125, -7., 'Sgr Le', charsize = 0.7, color = 255

xyouts, 91, -69., 'LMC', charsize = 0.7, color = 255
xyouts, 27, -75., 'SMC', charsize = 0.7, color = 255

loadct, 0

print, ntar, ntar2, ntar3

area = float(tarx[1]-tarx[0])*(tary[1]-tary[0])
print, float(ntar)/area, float(ntar2)/area, float(ntar3)/area
print, !pi*float(ntar)/area, !pi*float(ntar2)/area, !pi*float(ntar3)/area

!p.multi = 0

device, /close
set_plot,'x'

;GOTO, finish

use = isection([wg, out], 2, n = nuse)
dist = 10^(0.2*(dmag+5))*1e-3

dg = {ra:(d.ra)[0], dec:(d.dec)[0], ebv:(d.ebv)[0], pmra:(d.pmra)[0], pmdec:(d.pmdec)[0], pmra_error:(d.pmra_error)[0], pmdec_error:(d.pmdec_error)[0], phot_g_mean_mag:(d.PHOT_G_MEAN_MAG)[0], bp_rp:(d.bp_rp)[0], parallax:(d.PARALLAX)[0], parallax_error:(d.parallax_error)[0], dist:dist[0], col:col[0], mg:mg[0], rpm:rpm[0], source_id:(d.source_id)[0]}

dg = replicate(dg, nuse)

dg.ra = (d.ra)[use]
dg.dec = (d.dec)[use]
dg.ebv = (d.ebv)[use]
dg.pmra = (d.pmra)[use]
dg.pmdec = (d.pmdec)[use]
dg.pmra_error = (d.pmra_error)[use]
dg.pmdec_error = (d.pmdec_error)[use]
dg.phot_g_mean_mag = (d.phot_g_mean_mag)[use]
dg.bp_rp = (d.bp_rp)[use]
dg.parallax = (d.parallax)[use]
dg.parallax_error = (d.parallax_error)[use]
dg.dist = dist[use]
dg.col = col[use]
dg.mg = mg[use]
dg.rpm = rpm[use]

dg.source_id = (d.source_id)[use]

path = '~/Documents/Work/lists/rg_gaia2mass.fits'
mwrfits, dg, path, /create


;; snapshot = TVRD(True = 1)
;; write_png, '~/work/plots/gaia_edr3_lmc/pisces_giants.png', snapshot

finish:

END 
