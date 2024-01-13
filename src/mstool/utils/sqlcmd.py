# in nonbonded, param is the atom type in integer

sql_create = """

CREATE TABLE global_cell (id integer primary key, x float, y float, z float);

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
  segname text not null,
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

CREATE TABLE bond (p0 integer, p1 integer, 'order' integer);








CREATE TABLE stretch_harm_term (p0 integer, p1 integer, 'constrained' integer, param integer not null);
CREATE TABLE stretch_harm_param ('type' text, 'r0' float, 'fc' float, 'memo' text, id integer primary key);

CREATE TABLE angle_harm_term (p0 integer, p1 integer, p2 integer, 'constrained' integer, param integer not null);
CREATE TABLE angle_harm_param ('type' text, 'theta0' float, 'fc' float, 'memo' text, id integer primary key);

CREATE TABLE angle_harmcos_term (p0 integer, p1 integer, p2 integer, param integer not null);
CREATE TABLE angle_harmcos_param ('type' text, 'cos_theta0' float, 'fc' float, 'memo' text, id integer primary key);

CREATE TABLE dihedral_trig_term (p0 integer, p1 integer, p2 integer, p3 integer, param integer not null);
CREATE TABLE dihedral_trig_param ('type' text, 'phi0' float, 'fc0' float, 'fc1' float, 'fc2' float, 'fc3' float, 'fc4' float, 'fc5' float, 'fc6' float, 'memo' text, id integer primary key);

CREATE TABLE improper_harm_param ('type' text, 'phi0' float, 'fc' float, 'memo' text, id integer primary key);
CREATE TABLE improper_harm_term (p0 integer, p1 integer, p2 integer, p3 integer, param integer not null);

CREATE TABLE pair_12_6_es_param ('aij' float, 'bij' float, 'qij' float, 'type' text, 'memo' text, id integer primary key);
CREATE TABLE pair_12_6_es_term (p0 integer, p1 integer, param integer not null);

CREATE TABLE posre_harm_param ('fcx' float, 'fcy' float, 'fcz' float, id integer primary key);
CREATE TABLE posre_harm_term (p0 integer, 'x0' float, 'y0' float, 'z0' float, param integer not null);


CREATE TABLE bond_term (name text);
INSERT INTO "bond_term" VALUES('stretch_harm');
INSERT INTO "bond_term" VALUES('angle_harm');
INSERT INTO "bond_term" VALUES('angle_harmcos');
INSERT INTO "bond_term" VALUES('dihedral_trig');
INSERT INTO "bond_term" VALUES('improper_harm');
INSERT INTO "bond_term" VALUES('pair_12_6_es');
INSERT INTO "bond_term" VALUES('posre_harm');


CREATE VIEW stretch_harm as
  select p0, p1,
"type", "r0", "fc", "memo", "constrained"  from stretch_harm_param
  join stretch_harm_term
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
CREATE VIEW posre_harm as
  select p0,
"fcx", "fcy", "fcz", "x0", "y0", "z0"  from posre_harm_param
  join posre_harm_term
  on param=id;









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







CREATE TABLE nonbonded_info (vdw_funct text, vdw_rule text, es_funct text);
INSERT INTO  nonbonded_info VALUES('vdw_12_6','arithmetic/geometric','');
/*CREATE TABLE nonbonded_param ('type' text, 'sigma' float, 'epsilon' float, 'nbfix_identifier' text, 'memo' text, id integer primary key);*/
CREATE TABLE nonbonded_param ('type' text, 'sigma' float, 'epsilon' float, id integer primary key);
/*CREATE TABLE nonbonded_combined_param ('param1' integer, 'param2' integer, 'type' text, 'sigma' float, 'epsilon' float, 'nbfix_identifier' text, 'memo' text);*/
CREATE TABLE nonbonded_combined_param ('param1' integer, 'param2' integer, 'type' text, 'sigma' float, 'epsilon' float);




CREATE TABLE exclusion (p0 integer, p1 integer);
CREATE TABLE provenance (id integer primary key, version text, timestamp text, user text, workdir text, cmdline text, executable text);
"""



sql_insert_particle = """
INSERT INTO particle
    (id, anum, name, resname, chain, resid, mass, charge, x, y, z, vx, vy, vz, segname, insertion, msys_ct, nbtype, type, bfactor)
    VALUES
    ('{:d}', '{:d}', '{:s}', '{:s}', '{:s}', '{:d}', '{:f}', '{:f}', '{:f}', '{:f}', '{:f}', '{:f}', '{:f}', '{:f}', '{:s}', '{:s}', '{:d}', '{:d}', '{:s}', '{:f}');
"""

sql_insert_particle = """
INSERT INTO particle
    (id, anum, name, resname, chain, resid, mass, charge, x, y, z, vx, vy, vz, segname, insertion, msys_ct, nbtype, type, bfactor)
    VALUES
    (?,  ?,    ?,    ?,       ?,     ?,     ?,    ?,      ?, ?, ?, ?,  ?,  ?,  ?,       ?,         ?,       ?,      ?,    ?);
"""

sql_insert_exclusion = "INSERT INTO exclusion VALUES('{:d}','{:d}');"
sql_insert_bond      = "INSERT INTO bond      VALUES('{:d}', '{:d}', '{:d}');"

# fcx, fcy, fcz, posre index
sql_insert_posre_harm_param = "INSERT INTO posre_harm_param VALUES('{:f}','{:f}','{:f}',{:d});"
# atom index, x0, y0, z0, posre index
sql_insert_posre_harm_term  = "INSERT INTO posre_harm_term  VALUES('{:d}','{:f}','{:f}','{:f}','{:d}');"

sql_insert_stretch_harm_term = """
  INSERT INTO 'stretch_harm_term' 
  (p0, p1, 'constrained', param)
  VALUES
  ('{:d}', '{:d}', '{:d}', '{:d}');
"""

sql_insert_stretch_harm_param = """
  INSERT INTO 'stretch_harm_param'
  ('type', 'r0', 'fc', 'memo', id)
  VALUES
  ('{:s}', '{:f}', '{:f}', '', '{:d}');
"""

sql_insert_angle_harm_term = """
  INSERT INTO 'angle_harm_term'
  (p0, p1, p2, 'constrained', param)
  VALUES
  ('{:d}', '{:d}', '{:d}', '{:d}', '{:d}');
"""

sql_insert_angle_harm_param = """
  INSERT INTO 'angle_harm_param'
  ('type', 'theta0', 'fc', 'memo', id)
  VALUES
  ('{:s}', '{:f}', '{:f}', '', '{:d}');
"""

sql_insert_angle_harmcos_term = """
  INSERT INTO 'angle_harmcos_term'
  (p0, p1, p2, param)
  VALUES
  ('{:d}', '{:d}', '{:d}', '{:d}');
"""

sql_insert_angle_harmcos_param = """
  INSERT INTO 'angle_harmcos_param'
  ('type', 'cos_theta0', 'fc', 'memo', id)
  VALUES
  ('{:s}', '{:f}', '{:f}', '', '{:d}');
"""

sql_insert_dihedral_trig_term = """
 INSERT INTO 'dihedral_trig_term'
 (p0, p1, p2, p3, param)
 VALUES
 ('{:d}', '{:d}', '{:d}', '{:d}', '{:d}');
"""

sql_insert_dihedral_trig_param = """
  INSERT INTO 'dihedral_trig_param'
  ('type', 'phi0', 'fc0', 'fc1', 'fc2', 'fc3', 'fc4', 'fc5', 'fc6', 'memo', id)
  VALUES
  ('{:s}', '{:f}', '{:f}','{:f}','{:f}','{:f}','{:f}','{:f}','{:f}', '', '{:d}');
"""

sql_insert_vlc3_param = """
  INSERT INTO 'virtual_lc3_param' 
  ('c1', 'c2', id)
  VALUES
  ('{:f}', '{:f}', '{:d}');
"""

sql_insert_vlc3_term = """
  INSERT INTO 'virtual_lc3_term'
  ('p0', 'p1', 'p2', 'p3', param)
  VALUES
  ('{:d}', '{:d}', '{:d}', '{:d}', '{:d}');
"""

sql_insert_vout3_param = """
  INSERT INTO 'virtual_out3_param' 
  ('c1', 'c2', 'c3', id)
  VALUES
  ('{:f}', '{:f}', '{:f}', '{:d}');
"""

sql_insert_vout3_term = """
  INSERT INTO 'virtual_out3_term'
  ('p0', 'p1', 'p2', 'p3', param)
  VALUES
  ('{:d}', '{:d}', '{:d}', '{:d}', '{:d}');
"""

sql_insert_improper_harm_term = """
  INSERT INTO 'improper_harm_term'
  (p0, p1, p2, p3, param)
  VALUES
  ('{:d}', '{:d}', '{:d}', '{:d}', '{:d}')
"""

sql_insert_improper_harm_param = """
  INSERT INTO 'improper_harm_param'
  ('type', 'phi0', 'fc', 'memo', id)
  VALUES
  ('{:s}', '{:f}', '{:f}', '', '{:d}');
"""

sql_insert_cell = "INSERT INTO 'global_cell' VALUES({:d}, {:f}, {:f}, {:f});"


sql_insert_nonbonded = """
INSERT INTO nonbonded_param
    (id, epsilon, sigma, type)
    VALUES
    ('{:d}', '{:f}', '{:f}', '{:s}')
"""

### param is the atom type in integer
sql_insert_nonbonded_combined = """
INSERT INTO nonbonded_combined_param
    (param1, param2, epsilon, sigma, type)
    VALUES
    ('{:d}', '{:d}', '{:f}', '{:f}', '{:s}')
"""
