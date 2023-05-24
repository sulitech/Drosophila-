-- MySQL dump 10.13  Distrib 5.7.37, for Win64 (x86_64)
--
-- Host: localhost    Database: bioinfo
-- ------------------------------------------------------
-- Server version	5.7.37-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `tissue`
--

DROP TABLE IF EXISTS `tissue`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tissue` (
  `tissueid` int(12) NOT NULL COMMENT 'tissueid',
  `tissuename` varchar(200) NOT NULL COMMENT 'tissue name',
  `tname` varchar(200) NOT NULL,
  `genecount` int(12) NOT NULL DEFAULT '1',
  `cellcount` int(12) NOT NULL DEFAULT '0',
  `index_1` int(12) DEFAULT NULL,
  `fullpath` varchar(300) DEFAULT NULL,
  `dispord` int(12) NOT NULL,
  PRIMARY KEY (`tissueid`),
  KEY `index_partname` (`tissuename`),
  KEY `index_pname` (`tname`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `tissue`
--

LOCK TABLES `tissue` WRITE;
/*!40000 ALTER TABLE `tissue` DISABLE KEYS */;
INSERT INTO `tissue` VALUES (1,'Testis','Testis',15833,44621,NULL,'D:/project/python/bioinfo/data/loom/r_fca_biohub_testis_10x.loom',15),(2,'Antenna','Antenna',13203,37254,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_antenna_10x.loom',3),(3,'Body','Body',15267,96926,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_body_10x.loom',2),(4,'Body_wall','Body wall',12924,16551,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_body_wall_10x.loom',4),(5,'Fat_body','Fat body',15524,26926,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_fat_body_10x.loom',5),(6,'Gut','Gut',13407,11788,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_gut_10x.loom',6),(7,'Haltere','Haltere',12158,6527,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_haltere_10x.loom',7),(8,'Head','Head',13056,100527,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_head_10x.loom',1),(9,'Heart','Heart',12674,10686,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_heart_10x.loom',8),(10,'Leg','Leg',11164,14197,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_leg_10x.loom',9),(11,'Male_reproductive_glands','Male reproductive glands',13668,13143,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_male_reproductive_glands_10x.loom',10),(12,'Malpighian_tubule','Malpighian tubule',12843,13774,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_malpighian_tubule_10x.loom',11),(13,'Oenocyte','Oenocyte',14644,14410,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_oenocyte_10x.loom',12),(14,'Ovary','Ovary',12460,31401,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_ovary_10x.loom',13),(15,'proboscis_and_maxillary_palps','Proboscis',12802,26301,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_proboscis_and_maxillary_palps_10x.loom',14),(16,'Trachea','Trachea',14649,26906,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_trachea_10x.loom',16),(17,'Wing','Wing',13411,15889,NULL,'D:/project/python/bioinfo/data/loom/s_fca_biohub_wing_10x.loom',17);
/*!40000 ALTER TABLE `tissue` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2023-05-01 10:35:27
