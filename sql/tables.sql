drop database depths_exome_5bp;
create database depths_exome_5bp;
use depths_exome_5bp;


CREATE TABLE sample (

  sid                 INT NOT NULL AUTO_INCREMENT PRIMARY KEY,

  name                VARCHAR(80) NOT NULL,
  total_reads	      INT,
  mapped_reads	      INT,
  duplicate_reads     INT,
  mean_isize	      float
  gender              EMUM ('F', 'M', 'U' ) DEFAULT 'U',

);


CREATE TABLE region (

  rid                 INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  chr                 VARCHAR(8) NOT NULL ,
  start               INT NOT NULL,
  end                 INT NOT NULL,

  KEY pos_idx  (chr, start, end)
);

CREATE TABLE transcript (
 rid INT NOT NULL,

 exon_name varchar(300),
 refseq    varchar(300),
 CCDS      varchar(300),

 KEY ref_idx (refseq),
 KEY name_idx (exon_name),

 KEY rid_idx (rid)
);
 

CREATE TABLE coverage (

  sid                INT NOT NULL,
  rid                 INT NOT NULL,

  min		      FLOAT,  
  mean		      FLOAT,  
  max		      FLOAT,  
  missing	      VARCHAR(300),
  1to5		      VARCHAR(300),
  6to9		      VARCHAR(300),
  10to19	      VARCHAR(300),

  KEY sid_idx  (sid),
  KEY rid_idx  (rid)
);



