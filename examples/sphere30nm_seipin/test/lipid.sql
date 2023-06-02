BEGIN TRANSACTION;
CREATE TABLE angle_harm_param ('type' text, 'theta0' float, 'fc' float, 'memo' text, id integer primary key);
CREATE TABLE angle_harm_term (p0 integer, p1 integer, p2 integer, 'constrained' integer, param integer not null);
CREATE TABLE angle_harmcos_param ('type' text, 'cos_theta0' float, 'fc' float, 'memo' text, id integer primary key);
CREATE TABLE angle_harmcos_term (p0 integer, p1 integer, p2 integer, param integer not null);
CREATE TABLE bond (p0 integer, p1 integer, 'order' integer);
CREATE TABLE bond_term (name text);
INSERT INTO "bond_term" VALUES('stretch_harm');
INSERT INTO "bond_term" VALUES('posre_harm');
INSERT INTO "bond_term" VALUES('angle_harm');
INSERT INTO "bond_term" VALUES('angle_harmcos');
INSERT INTO "bond_term" VALUES('dihedral_trig');
INSERT INTO "bond_term" VALUES('improper_harm');
INSERT INTO "bond_term" VALUES('pair_12_6_es');
CREATE TABLE dihedral_trig_param ('type' text, 'phi0' float, 'fc0' float, 'fc1' float, 'fc2' float, 'fc3' float, 'fc4' float, 'fc5' float, 'fc6' float, 'memo' text, id integer primary key);
CREATE TABLE dihedral_trig_term (p0 integer, p1 integer, p2 integer, p3 integer, param integer not null);
CREATE TABLE exclusion (p0 integer, p1 integer);
CREATE TABLE global_cell (id integer primary key, x float, y float, z float);
INSERT INTO "global_cell" VALUES(1,0.0,0.0,0.0);
INSERT INTO "global_cell" VALUES(2,0.0,0.0,0.0);
INSERT INTO "global_cell" VALUES(3,0.0,0.0,0.0);
CREATE TABLE improper_harm_param ('type' text, 'phi0' float, 'fc' float, 'memo' text, id integer primary key);
CREATE TABLE improper_harm_term (p0 integer, p1 integer, p2 integer, p3 integer, param integer not null);
CREATE TABLE nonbonded_combined_param ('param1' integer, 'param2' integer, 'type' text, 'sigma' float, 'epsilon' float);
CREATE TABLE nonbonded_info (vdw_funct text, vdw_rule text, es_funct text);
INSERT INTO "nonbonded_info" VALUES('vdw_12_6','arithmetic/geometric','');
CREATE TABLE nonbonded_param ('type' text, 'sigma' float, 'epsilon' float, id integer primary key);
CREATE TABLE pair_12_6_es_param ('aij' float, 'bij' float, 'qij' float, 'type' text, 'memo' text, id integer primary key);
CREATE TABLE pair_12_6_es_term (p0 integer, p1 integer, param integer not null);
CREATE TABLE particle (
  id integer primary key,
  anum integer,
  name text not null,
  x float,
  y float,
  z float,
  vx float,
  vy float,
  vz float,
  resname text not null,
  resid integer,
  chain text not null,
  segid text not null,
  mass float,
  charge float,
  formal_charge integer,
  insertion text not null,
  msys_ct integer not null,
  'm_grow_name' text,
  'm_mmod_type' integer,
  nbtype integer not null,
  type text not null,
  bfactor float
);
INSERT INTO "particle" VALUES(0,7,'N',-32.066,33.682,23.148,0.0,0.0,0.0,'POPC',101,'X','X',14.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(1,6,'C12',-32.732,32.612,22.283,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(2,1,'H12A',-32.98,31.782,22.93,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(3,1,'H12B',-33.666,33.027,21.927,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(4,6,'C13',-30.845,33.135,23.852,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(5,1,'H13A',-31.127,32.288,24.456,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(6,1,'H13B',-30.415,33.904,24.474,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(7,1,'H13C',-30.119,32.822,23.116,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(8,6,'C14',-31.683,34.885,22.315,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(9,1,'H14A',-31.238,35.638,22.949,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(10,1,'H14B',-32.567,35.285,21.838,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(11,1,'H14C',-30.967,34.583,21.556,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(12,6,'C15',-33.061,34.126,24.195,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(13,1,'H15A',-32.62,34.891,24.815,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(14,1,'H15B',-33.338,33.279,24.804,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(15,1,'H15C',-33.943,34.516,23.705,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(16,6,'C11',-31.985,32.038,21.042,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(17,1,'H11A',-32.649,31.289,20.557,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(18,1,'H11B',-31.816,32.865,20.321,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(19,15,'P',-29.446,31.725,20.587,0.0,0.0,0.0,'POPC',101,'X','X',31.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(20,8,'O13',-29.735,31.735,19.134,0.0,0.0,0.0,'POPC',101,'X','X',16.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(21,8,'O14',-28.375,30.821,21.063,0.0,0.0,0.0,'POPC',101,'X','X',16.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(22,8,'O12',-30.77,31.396,21.398,0.0,0.0,0.0,'POPC',101,'X','X',16.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(23,8,'O11',-29.133,33.208,21.069,0.0,0.0,0.0,'POPC',101,'X','X',16.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(24,6,'C1',-27.914,33.769,20.596,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(25,1,'HA',-27.068,33.438,21.242,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(26,1,'HB',-27.701,33.432,19.555,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(27,6,'C2',-27.983,35.327,20.604,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(28,1,'HS',-27.851,35.63,21.669,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(29,8,'O21',-29.26,35.849,20.181,0.0,0.0,0.0,'POPC',101,'X','X',16.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(30,6,'C21',-29.83,35.526,19.014,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(31,8,'O22',-29.393,34.763,18.17,0.0,0.0,0.0,'POPC',101,'X','X',16.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(32,6,'C22',-31.154,36.308,18.898,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(33,1,'H2R',-31.051,37.291,19.408,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(34,1,'H2S',-31.946,35.725,19.421,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(35,6,'C3',-26.82,35.968,19.795,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(36,1,'HX',-27.015,35.802,18.716,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(37,1,'HY',-25.874,35.444,20.039,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(38,8,'O31',-26.656,37.352,20.138,0.0,0.0,0.0,'POPC',101,'X','X',16.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(39,6,'C31',-26.45,38.361,19.262,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(40,8,'O32',-26.003,39.442,19.615,0.0,0.0,0.0,'POPC',101,'X','X',16.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(41,6,'C32',-26.864,38.064,17.792,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(42,1,'H2X',-26.16,37.286,17.423,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(43,1,'H2Y',-27.897,37.649,17.784,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(44,6,'C23',-31.591,36.543,17.431,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(45,1,'H3R',-32.569,37.076,17.469,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(46,1,'H3S',-31.758,35.569,16.927,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(47,6,'C24',-30.597,37.386,16.604,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(48,1,'H4R',-29.685,36.769,16.443,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(49,1,'H4S',-30.293,38.278,17.194,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(50,6,'C25',-31.086,37.847,15.212,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(51,1,'H5R',-30.168,37.999,14.595,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(52,1,'H5S',-31.682,37.056,14.709,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(53,6,'C26',-31.831,39.201,15.193,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(54,1,'H6R',-31.397,39.853,15.977,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(55,1,'H6S',-31.618,39.694,14.218,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(56,6,'C27',-33.363,39.13,15.391,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(57,1,'H7R',-33.646,39.949,16.091,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(58,1,'H7S',-33.651,38.174,15.883,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(59,6,'C28',-34.214,39.326,14.111,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(60,1,'H8R',-34.018,40.357,13.723,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(61,1,'H8S',-35.279,39.305,14.43,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(62,6,'C29',-33.999,38.289,13.031,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(63,1,'H91',-33.668,37.302,13.399,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(64,6,'C210',-34.176,38.44,11.704,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(65,1,'H101',-33.99,37.57,11.059,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(66,6,'C211',-34.626,39.687,10.977,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(67,1,'H11R',-34.502,40.589,11.616,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(68,1,'H11S',-35.71,39.613,10.74,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(69,6,'C212',-33.861,39.964,9.661,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(70,1,'H12R',-32.767,39.939,9.851,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(71,1,'H12S',-34.119,40.997,9.339,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(72,6,'C213',-34.195,39.005,8.5,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(73,1,'H13R',-35.298,38.986,8.353,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(74,1,'H13S',-33.867,37.975,8.776,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(75,6,'C214',-33.531,39.409,7.169,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(76,1,'H14R',-33.843,40.453,6.935,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(77,1,'H14S',-32.429,39.409,7.296,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(78,6,'C215',-33.902,38.528,5.96,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(79,1,'H15R',-33.4,38.964,5.063,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(80,1,'H15S',-34.997,38.588,5.797,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(81,6,'C216',-33.495,37.044,6.084,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(82,1,'H16R',-32.403,36.983,6.276,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(83,1,'H16S',-34.023,36.602,6.961,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(84,6,'C217',-33.826,36.176,4.856,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(85,1,'H17R',-34.924,36.227,4.675,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(86,1,'H17S',-33.589,35.125,5.112,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(87,6,'C218',-33.093,36.527,3.556,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(88,1,'H18R',-31.989,36.535,3.715,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(89,1,'H18S',-33.404,37.528,3.184,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(90,1,'H18T',-33.325,35.773,2.772,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(91,6,'C33',-26.804,39.328,16.807,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(92,1,'H3X',-27.53,39.754,17.566,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(93,1,'H3Y',-26.849,38.482,16.055,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(94,6,'C34',-27.896,39.954,15.824,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(95,1,'H4X',-28.102,41.0,16.141,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(96,1,'H4Y',-28.839,39.385,15.91,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(97,6,'C35',-27.491,40.023,14.329,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(98,1,'H5X',-26.492,40.503,14.275,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(99,1,'H5Y',-27.407,38.988,13.929,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(100,6,'C36',-28.44,40.84,13.438,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(101,1,'H6X',-29.473,40.43,13.495,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(102,1,'H6Y',-28.472,41.872,13.856,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(103,6,'C37',-27.986,40.955,11.966,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(104,1,'H7X',-26.883,40.918,11.873,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(105,1,'H7Y',-28.295,41.974,11.63,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(106,6,'C38',-28.657,39.986,10.972,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(107,1,'H8X',-29.724,39.84,11.255,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(108,1,'H8Y',-28.64,40.495,9.983,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(109,6,'C39',-27.985,38.613,10.783,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(110,1,'H9X',-27.88,38.13,11.774,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(111,1,'H9Y',-26.962,38.754,10.378,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(112,6,'C310',-28.77,37.658,9.852,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(113,1,'H10X',-29.765,37.468,10.317,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(114,1,'H10Y',-28.224,36.688,9.808,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(115,6,'C311',-28.999,38.174,8.412,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(116,1,'H11X',-28.022,38.286,7.898,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(117,1,'H11Y',-29.468,39.184,8.474,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(118,6,'C312',-29.947,37.318,7.531,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(119,1,'H12X',-30.912,37.181,8.06,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(120,1,'H12Y',-30.146,37.948,6.633,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(121,6,'C313',-29.38,35.899,7.068,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(122,1,'H13X',-30.144,36.06,6.249,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(123,1,'H13Y',-29.415,35.728,8.184,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(124,6,'C314',-30.07,34.459,7.087,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(125,1,'H14X',-29.635,33.846,7.91,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(126,1,'H14Y',-31.156,34.557,7.285,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(127,6,'C315',-29.844,33.6,5.822,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(128,1,'H15X',-28.753,33.604,5.617,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(129,1,'H15Y',-30.13,32.551,6.068,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(130,6,'C316',-30.587,34.007,4.542,0.0,0.0,0.0,'POPC',101,'X','X',12.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(131,1,'H16X',-31.683,33.978,4.708,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(132,1,'H16Y',-30.292,35.031,4.231,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
INSERT INTO "particle" VALUES(133,1,'H16Z',-30.336,33.307,3.718,0.0,0.0,0.0,'POPC',101,'X','X',1.0,0.0,NULL,'',0,NULL,NULL,0,'tbd',0.0);
CREATE TABLE posre_harm_param ('fcx' float, 'fcy' float, 'fcz' float, id integer primary key);
CREATE TABLE posre_harm_term (p0 integer, 'x0' float, 'y0' float, 'z0' float, param integer not null);
CREATE TABLE provenance (id integer primary key, version text, timestamp text, user text, workdir text, cmdline text, executable text);
CREATE TABLE stretch_harm_param ('type' text, 'r0' float, 'fc' float, 'memo' text, id integer primary key);
CREATE TABLE stretch_harm_term (p0 integer, p1 integer, 'constrained' integer, param integer not null);
CREATE TABLE virtual_lc3_param ('c1' float, 'c2' float, id integer primary key);
CREATE TABLE virtual_lc3_term (p0 integer, p1 integer, p2 integer, p3 integer, param integer not null);
CREATE TABLE virtual_out3_param ('c1' float, 'c2' float, 'c3' float, id integer primary key);
CREATE TABLE virtual_out3_term (p0 integer, p1 integer, p2 integer, p3 integer, param integer not null);
CREATE TABLE virtual_term (name text);
INSERT INTO "virtual_term" VALUES('virtual_lc3');
INSERT INTO "virtual_term" VALUES('virtual_out3');
CREATE VIEW virtual_lc3 as
  select p0, p1, p2, p3,
"c1", "c2"  from virtual_lc3_param
  join virtual_lc3_term
  on param=id;
CREATE VIEW virtual_out3 as
  select p0, p1, p2, p3,
"c1", "c2", "c3"  from virtual_out3_param
  join virtual_out3_term
  on param=id;
CREATE VIEW angle_harm as
  select p0, p1, p2,
"type", "theta0", "fc", "memo", "constrained"  from angle_harm_param
  join angle_harm_term
  on param=id;
CREATE VIEW angle_harmcos as
  select p0, p1, p2,
"type", "cos_theta0", "fc", "memo" from angle_harmcos_param
  join angle_harmcos_term
  on param=id;
CREATE VIEW dihedral_trig as
  select p0, p1, p2, p3,
"type", "phi0", "fc0", "fc1", "fc2", "fc3", "fc4", "fc5", "fc6", "memo"  from dihedral_trig_param
  join dihedral_trig_term
  on param=id;
CREATE VIEW improper_harm as
  select p0, p1, p2, p3,
"type", "phi0", "fc", "memo"  from improper_harm_param
  join improper_harm_term
  on param=id;
CREATE VIEW pair_12_6_es as
  select p0, p1,
"aij", "bij", "qij", "type", "memo"  from pair_12_6_es_param
  join pair_12_6_es_term
  on param=id;
CREATE VIEW stretch_harm as
  select p0, p1,
"type", "r0", "fc", "memo", "constrained"  from stretch_harm_param
  join stretch_harm_term
  on param=id;
COMMIT;
