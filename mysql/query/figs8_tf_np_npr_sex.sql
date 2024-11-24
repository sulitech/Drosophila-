
---np

insert into prt_10_sex_exp (tname, motifname, tfgenename, tfcount, npgenename, tfauc, 
npmotifregoccur, npmotifregweight, npgexpress, nprtype, sex) 
select e.tname , np.motifname,  b.genename  tfgenename , count(np.tfgeneid) tfcount , 
c.genename  npgenename  , sum(np.motifauc) , sum(np.motifregoccur) , sum(np.motifregweight), sum(np.geneexpress),
 'np', np.sex
from 
(
		select b.tissueid , a.cellid , b.tfgeneid,  b.motifname ,  d.geneid , i.geneexpress, e.pid ,  b.motifauc ,
		c.reggeneid npgeneid , c.motifregoccur , c.motifregweight , a.sex
		from 
		(
			select cellid , sex from np_10_cell
		) a,
		motifregulonauc b ,
		motifregulon c ,
		nprtgene d ,
		pairprgene e,
		tissue_cell_gene i
		where a.cellid = b.cellid
		and b.tissueid=17
		and a.cellid = i.cellid
		and b.tissueid = c.tissueid
		and c.reggeneid = d.geneid
		and d.geneid = i.geneid
		and c.tfgeneid = b.tfgeneid
		and d.used=1
		and d.genetype='p'
		and b.tissueid= i.tissueid
		and d.geneid= e.geneid
) np ,
cell a ,
gene b , 
gene c ,
tissue e 
where np.cellid = a.cellid
and np.tfgeneid = b.geneid
and np.npgeneid = c.geneid
and np.tissueid = e.tissueid
group by e.tname , a.annotation, np.motifname, np.tfgeneid , np.npgeneid , np.sex
order by e.dispord , a.annotation , np.motifname, np.tfgeneid , np.npgeneid , np.sex

--npr

select count(*) from prt_10_exp

select * from prt_10_sex_exp 

insert into prt_10_sex_exp (tname, motifname, tfgenename, tfcount, nprgenename, tfauc, 
nprmotifregoccur, nprmotifregweight, nprgexpress, nprtype, sex) 
select e.tname , np.motifname,  b.genename  tfgenename , count(np.tfgeneid) tfcount , 
c.genename  nprgenename  ,  sum(np.motifauc) ,  sum(np.motifregoccur) , sum(np.motifregweight) 
,sum(np.geneexpress), 'npr', np.sex
from 
(
		select b.tissueid , a.cellid , b.tfgeneid,  b.motifname ,  d.geneid , i.geneexpress, e.pid ,  b.motifauc ,
		c.reggeneid npgeneid , c.motifregoccur , c.motifregweight , a.sex
		from 
		(
			select cellid , sex from npr_10_cell
		) a,
		motifregulonauc b ,
		motifregulon c ,
		nprtgene d ,
		pairprgene e,
		tissue_cell_gene i
		where a.cellid = b.cellid
		and b.tissueid=17
		and a.cellid = i.cellid
		and b.tissueid = c.tissueid
		and c.reggeneid = d.geneid
		and d.geneid = i.geneid
		and c.tfgeneid = b.tfgeneid
		and d.used=1
		and d.genetype='r'
		and b.tissueid= i.tissueid
		and d.geneid= e.geneid
) np ,
cell a ,
gene b , 
gene c ,
tissue e 
where np.cellid = a.cellid
and np.tfgeneid = b.geneid
and np.npgeneid = c.geneid
and np.tissueid = e.tissueid
group by e.tname , a.annotation, np.motifname, np.tfgeneid , np.npgeneid , np.sex
order by e.dispord , a.annotation , np.motifname, np.tfgeneid , np.npgeneid , np.sex


select count(b.aa) from 
(
select count(nprgenename) aa from prt_10_sex_exp
where nprtype='npr'
group by motifname, tfgenename, tfcount, npgenename
) b

select count(*) from 
(
select motifname, tfgenename, npgenename from prt_10_sex_exp
where nprtype='np'
group by motifname, tfgenename, npgenename
) b

select count(*) from 
(
select motifname, tfgenename, sum(tfcount), nprgenename , sex from prt_10_sex_exp
where nprtype='npr'
group by motifname, tfgenename, nprgenename , sex
) b

select count(*) from 
(
select motifname, tfgenename, sum(tfcount), npgenename , sex from prt_10_sex_exp
where nprtype='np' and sex='male'
group by motifname, tfgenename, npgenename , sex
) b


select * from prt_10_sex_exp
where tname='wing'
where sex='mix'


--all

    select replace(motifname,'motif','NP') , sum(tfcount), npgenename gname , sum(npmotifregoccur) gcount , nprtype ,sex
    from 
    (
		select a.* from prt_10_sex_exp a, tissue c
		where  a.nprtype='np' and a.sex='female'
		and a.tname = c.tname
		and not exists
		( 
		   select 1 from prt_10_exp b
			where nprtype='np_npr'
			and b.motifname=a.motifname
		) 
		order by c.dispord , a.motifname, npgenename
    ) aa
    group by aa.motifname , aa.npgenename , aa.nprtype
    
    	select replace(motifname,'motif','NPR') , sum(tfcount), nprgenename gname , sum(nprmotifregoccur) gcount , nprtype , sex
    from 
    (
		select a.* from prt_10_sex_exp a, tissue c
		where  a.nprtype='npr' and a.sex='female'
		and a.tname = c.tname
		and not exists
		( 
		    select 1 from prt_10_exp b
			where nprtype='np_npr'
			and b.motifname=a.motifname
		) 
		order by c.dispord , a.motifname, npgenename
    ) aa
    group by aa.motifname , aa.nprgenename , aa.nprtype

    
    
	select replace(motifname,'motif','NP') , sum(tfcount), npgenename gname , sum(npmotifregoccur) gcount , nprtype, sex
    from 
    (
		select a.* from prt_10_sex_exp a, tissue c
		where  a.nprtype='np' and a.sex='male'
		and a.tname = c.tname
		and not exists
		( 
		   select 1 from prt_10_exp b
			where nprtype='np_npr'
			and b.motifname=a.motifname
		) 
		order by c.dispord , a.motifname, npgenename
    ) aa
    group by aa.motifname , aa.npgenename , aa.nprtype
    
    

	select replace(motifname,'motif','NPR') , sum(tfcount), nprgenename gname , sum(nprmotifregoccur) gcount , nprtype ,sex
    from 
    (
		select a.* from prt_10_sex_exp a, tissue c
		where  a.nprtype='npr' and a.sex='male'
		and a.tname = c.tname
		and not exists
		( 
		   select 1 from prt_10_exp b
			where nprtype='np_npr'
			and b.motifname=a.motifname
		)        
		order by c.dispord , a.motifname, npgenename
    ) aa
    group by aa.motifname , aa.nprgenename , aa.nprtype 
    




select * from prt_10_sex_exp
--whole tissue
select a.tname , replace(a.motifname,'motif','NP') , sum(a.tfcount), a.npgenename gname , sum(a.tfauc), 
sum(a.npmotifregoccur) gcount ,  sum(a.npmotifregweight) ,  sum(a.npgexpress) ,  a.nprtype , a.sex
from prt_10_sex_exp a, tissue c
where  a.nprtype='np' and a.sex='female'
and a.tname = c.tname
group by  a.tname, a.motifname , a.npgenename , a.nprtype , a.sex
order by c.dispord , a.motifname, npgenename


select a.tname , replace(a.motifname,'motif','NPR') , sum(a.tfcount), a.nprgenename gname , sum(a.tfauc) ,  sum(a.nprmotifregoccur) gcount ,  
sum(a.nprmotifregweight) ,	sum(a.nprgexpress)  , a.nprtype , a.sex
from prt_10_sex_exp a, tissue c
where  a.nprtype='npr' and a.sex='female'
and a.tname = c.tname
group by  a.tname, a.motifname , a.nprgenename , a.nprtype , a.sex
order by c.dispord , a.motifname, a.nprgenename


select a.tname , replace(a.motifname,'motif','NP') , sum(a.tfcount), a.npgenename gname , sum(a.tfauc), 
sum(a.npmotifregoccur) gcount ,  sum(a.npmotifregweight) ,  sum(a.npgexpress) ,  a.nprtype , a.sex
from prt_10_sex_exp a, tissue c
where  a.nprtype='np' and a.sex='male'
and a.tname = c.tname
group by  a.tname, a.motifname , a.npgenename , a.nprtype , a.sex
order by c.dispord , a.motifname, npgenename


select a.tname ,replace(a.motifname,'motif','NPR') , sum(a.tfcount), a.nprgenename gname , sum(a.tfauc) ,  sum(a.nprmotifregoccur) gcount ,  
sum(a.nprmotifregweight) ,	sum(a.nprgexpress)  , a.nprtype , a.sex
from prt_10_sex_exp a, tissue c
where  a.nprtype='npr' and a.sex='male'
and a.tname = c.tname
group by  a.tname, a.motifname , a.nprgenename , a.nprtype , a.sex
order by c.dispord , a.motifname, a.nprgenename


--female만 존재

select replace(motifname,'motif','NP') , sum(tfcount), npgenename gname , sum(npmotifregoccur) gcount , nprtype, sex
from
(
select a.* from  prt_10_sex_exp a
where a.nprtype='np' and a.sex = 'female'
and not exists
( 
   select 1 from prt_10_sex_exp b
	where nprtype='np' and b.sex = 'male'
	and b.motifname=a.motifname
    and b.tname = a.tname 
    and b.npgenename = a.npgenename 
)   
) aa
group by aa.motifname , aa.npgenename , aa.nprtype 

select replace(motifname,'motif','NPR') , sum(tfcount), nprgenename gname , sum(nprmotifregoccur) gcount , nprtype, sex
from
(
select a.* from  prt_10_sex_exp a
where a.nprtype='npr' and a.sex = 'female'
and not exists
( 
   select 1 from prt_10_sex_exp b
	where nprtype='npr' and b.sex = 'male'
	and b.motifname=a.motifname
    and b.tname = a.tname 
    and b.nprgenename = a.nprgenename 
)         
) aa
group by aa.motifname , aa.nprgenename , aa.nprtype 
 

--male만 존재..
select replace(motifname,'motif','NP') , sum(tfcount), npgenename gname , sum(npmotifregoccur) gcount , nprtype, sex
from
(
select a.* from  prt_10_sex_exp a
where a.nprtype='np' and a.sex = 'male'
and not exists
( 
   select 1 from prt_10_sex_exp b
	where nprtype='np' and b.sex = 'female'
	and b.motifname=a.motifname
    and b.tname = a.tname 
    and b.npgenename = a.npgenename 
)   
) aa
group by aa.motifname , aa.npgenename , aa.nprtype 


select replace(motifname,'motif','NPR') , sum(tfcount), nprgenename gname , sum(nprmotifregoccur) gcount , nprtype, sex
from
(
select a.* from  prt_10_sex_exp a
where a.nprtype='npr' and a.sex = 'male'
and not exists
( 
   select 1 from prt_10_sex_exp b
	where nprtype='npr' and b.sex = 'female'
	and b.motifname=a.motifname
    and b.tname = a.tname 
    and b.nprgenename = a.nprgenename 
)          
) aa
group by aa.motifname , aa.nprgenename , aa.nprtype 
 

--whole tissue

select a.tname , replace(a.motifname,'motif','NP') , sum(a.tfcount), a.npgenename gname , sum(a.tfauc), 
sum(a.npmotifregoccur) gcount ,  sum(a.npmotifregweight) ,  sum(a.npgexpress) ,  a.nprtype , a.sex
from prt_10_sex_exp a, tissue c
where  a.nprtype='np' and a.sex='female'
and a.tname = c.tname
group by  a.tname, a.motifname , a.npgenename , a.nprtype , a.sex
order by c.dispord , a.motifname, npgenename


select a.tname , replace(a.motifname,'motif','NPR') , sum(a.tfcount), a.nprgenename gname , sum(a.tfauc) ,  sum(a.nprmotifregoccur) gcount ,  
sum(a.nprmotifregweight) ,	sum(a.nprgexpress)  , a.nprtype , a.sex
from prt_10_sex_exp a, tissue c
where  a.nprtype='npr' and a.sex='female'
and a.tname = c.tname
group by  a.tname, a.motifname , a.nprgenename , a.nprtype , a.sex
order by c.dispord , a.motifname, a.nprgenename





--whole tissue

select aa.tname ,  replace(aa.motifname,'motif','NP') , sum(aa.tfcount), aa.npgenename gname , sum(aa.tfauc) ,  
sum(aa.npmotifregoccur) gcount ,  sum(aa.npmotifregweight) ,	sum(aa.npgexpress)  , aa.nprtype , aa.sex
from
(
	select a.* from  prt_10_sex_exp a
	where a.nprtype='np' and a.sex = 'female'
	and not exists
	( 
	   select 1 from prt_10_sex_exp b
		where nprtype='np' and b.sex = 'male'
		and b.motifname=a.motifname
		and b.tname = a.tname 
		and b.npgenename = a.npgenename 
	)   
) aa
group by  aa.tname, aa.motifname , aa.npgenename , aa.nprtype , aa.sex
order by aa.tname, aa.motifname, aa.npgenename


select aa.tname , replace(aa.motifname,'motif','NPR') , sum(aa.tfcount), aa.nprgenename gname , sum(aa.nprmotifregoccur) gcount , sum(aa.tfauc) ,  
sum(aa.nprmotifregweight) ,	sum(aa.nprgexpress)  , aa.nprtype , aa.sex
from
(
select a.* from  prt_10_sex_exp a
where a.nprtype='npr' and a.sex = 'female'
and not exists
( 
   select 1 from prt_10_sex_exp b
	where nprtype='npr' and b.sex = 'male'
	and b.motifname=a.motifname
    and b.tname = a.tname 
    and b.nprgenename = a.nprgenename 
)         
) aa
group by aa.tname, aa.motifname , aa.nprgenename , aa.nprtype 
order by aa.tname , aa.motifname, aa.nprgenename
 

--male만 존재..


select aa.tname ,  replace(aa.motifname,'motif','NP') , sum(aa.tfcount), aa.npgenename gname , sum(aa.tfauc) ,  
sum(aa.npmotifregoccur) gcount ,  sum(aa.npmotifregweight) ,	sum(aa.npgexpress)  , aa.nprtype , aa.sex
from
(
select a.* from  prt_10_sex_exp a
where a.nprtype='np' and a.sex = 'male'
and not exists
( 
   select 1 from prt_10_sex_exp b
	where nprtype='np' and b.sex = 'female'
	and b.motifname=a.motifname
    and b.tname = a.tname 
    and b.npgenename = a.npgenename 
)   
) aa
group by  aa.tname, aa.motifname , aa.npgenename , aa.nprtype , aa.sex
order by aa.tname, aa.motifname, aa.npgenename


select aa.tname , replace(aa.motifname,'motif','NPR') , sum(aa.tfcount), aa.nprgenename gname , sum(aa.nprmotifregoccur) gcount , sum(aa.tfauc) ,  
sum(aa.nprmotifregweight) ,	sum(aa.nprgexpress)  , aa.nprtype , aa.sex
from
(
select a.* from  prt_10_sex_exp a
where a.nprtype='npr' and a.sex = 'male'
and not exists
( 
   select 1 from prt_10_sex_exp b
	where nprtype='npr' and b.sex = 'female'
	and b.motifname=a.motifname
    and b.tname = a.tname 
    and b.nprgenename = a.nprgenename 
)          
) aa
group by aa.tname, aa.motifname , aa.nprgenename , aa.nprtype 
order by aa.tname , aa.motifname, aa.nprgenename


ALTER TABLE ntr_8_tf_expr ADD COLUMN annotaion varchar(300) ;
ALTER TABLE npr_10_cell ADD COLUMN age varchar(10) ;

select count(*) from cell 
where age=2
1-44621,5 423137,3 26926,2 13143
where

select * from prt_10_sex_exp
where sex='male'
and tname='head'
and nprgenename like 'sNPF%'

select * from prt_10_sex_exp a
where a.sex='female'
and a.tname='head'
and a.npgenename ='sNPF'
and not exists
(
	select 1 from prt_10_sex_exp b
    where b.sex='female'
    and b.tname='head'
    and b.npgenename ='sNPF'
    and a.motifname = b.motifname
)


select distinct annotation from cell a
where a.sex='male'
and not exists
( select 1 from cell b
where b.sex='female'
and a.annotation = b.annotation

select * from 
	(
	select distinct annotation  from cell 
	where sex='male'
    and annotation not in 
    (
    	select distinct annotation  from cell 
	   where sex='female'
    )
	
select distinct a.annotation  from cell a,
tissue_cell b 
where a.cellid = b.cellid
and b.tissueid = 3
/*
a.annotation
and b.tissueid = 11*/
and a.annotation in
(
	select distinct annotation  from cell 
	where sex='female'
    and annotation not in 
    (
    	select distinct annotation  from cell 
	   where sex='male'
    )
)

