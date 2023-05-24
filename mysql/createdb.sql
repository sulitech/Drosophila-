--create database
CREATE DATABASE drosphila DEFAULT CHARACTER SET utf8 COLLATE utf8_general_ci;

-- create user root/admin
create user 'root'@'localhost' identified by 'admin';

--data import to mysql

mysql -uroot -padmin --default-character-set=utf8 drosphila < gene.sql
mysql -uroot -padmin --default-character-set=utf8 drosphila < cell.sql
mysql -uroot -padmin --default-character-set=utf8 drosphila < motifregulon.sql
mysql -uroot -padmin --default-character-set=utf8 drosphila < motifregulonauc.sql
mysql -uroot -padmin --default-character-set=utf8 drosphila < tissue.sql
mysql -uroot -padmin --default-character-set=utf8 drosphila < tissue_cell.sql
mysql -uroot -padmin --default-character-set=utf8 drosphila < tissue_cell_gene.sql
mysql -uroot -padmin --default-character-set=utf8 drosphila < tissue_gene.sql
