# SQL query to be pasted into the CAS entry box at
#   http://skyserver.sdss.org/dr12/en/tools/search/sql.aspx
#
# Based on the sample photo-z and LRG queries at 
#   http://skyserver.sdss.org/dr12/en/help/docs/realquery.aspx#photoz
#   http://skyserver.sdss.org/dr12/en/help/docs/realquery.aspx#lrg

SELECT
 ra, 
 dec,
 pz.z,
 dered_u,
 dered_g,
 dered_r,
 dered_i,
 dered_z
 FROM Galaxy AS g
 JOIN Photoz pz ON g.objid = pz.objid
 WHERE ( ( g.flags & (dbo.fPhotoFlags('BINNED1') 
 | dbo.fPhotoFlags('BINNED2') 
 | dbo.fPhotoFlags('BINNED4')) ) > 0 
 and ( g.flags & (dbo.fPhotoFlags('BLENDED') 
 | dbo.fPhotoFlags('NODEBLEND') 
 | dbo.fPhotoFlags('CHILD')) ) != dbo.fPhotoFlags('BLENDED') 
 and ( g.flags & (dbo.fPhotoFlags('EDGE') 
 | dbo.fPhotoFlags('SATURATED')) ) = 0 
 and g.petroMag_i > 17.5 
 and (g.petroMag_r > 15.5 or g.petroR50_r > 2) 
 and (g.petroMag_r > 0 and g.g > 0 and g.r > 0 and g.i > 0) 
 and ( (g.petroMag_r-g.extinction_r) < 19.2 
 and (g.petroMag_r - g.extinction_r < 
 (13.1 + (7/3) * (g.dered_g - g.dered_r) + 4 * (g.dered_r - g.dered_i) 
 - 4 * 0.18) ) 
 and ( (g.dered_r - g.dered_i - (g.dered_g - g.dered_r)/4 - 0.18) < 0.2) 
 and ( (g.dered_r - g.dered_i - (g.dered_g - g.dered_r)/4 - 0.18) > -0.2) 
 and ( (g.petroMag_r - g.extinction_r + 
 2.5 * LOG10(2 * 3.1415 * g.petroR50_r * g.petroR50_r)) < 24.2) ) 
 or ( (g.petroMag_r - g.extinction_r < 19.5) 
 and ( (g.dered_r - g.dered_i - (g.dered_g - g.dered_r)/4 - 0.18) > (0.45 - 4 * 
 (g.dered_g - g.dered_r)) ) 
 and ( (g.dered_g - g.dered_r) > (1.35 + 0.25 * (g.dered_r - g.dered_i)) ) ) 
 and ( (g.petroMag_r - g.extinction_r + 
 2.5 * LOG10(2 * 3.1415 * g.petroR50_r * g.petroR50_r) ) < 23.3 ) 
 and (pz.z BETWEEN 0.0 and 0.9)
 and (pz.photoErrorClass=1)
 and (pz.nnCount>95)
 and (pz.zErr BETWEEN 0 and 0.03 )
 and rand() <= 0.0001 )
