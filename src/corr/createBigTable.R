### get table ready for correlation

source('../bigVariateTable_join.R')

con <- dbCon()
metaTable <- dbGetQuery(con, sprintf('select * from nhanes.ewas_var_tab_combined where analyzable = 1'))
metaTable <- subset(metaTable, 
	var != 'bmi' & 
	var_desc_ewas != 'occupation' & 
	var_desc_ewas != 'food and supplement consumption' & 
	var_desc_ewas != 'food recall')


bt_a <- bigTable(con, 'a')
bt_b <- bigTable(con, 'b')
bt_c <- bigTable(con, 'c')
bt_d <- bigTable(con, 'd')

ocq_a <- dbGetQuery(con, 'select SEQN, occupation from nhanes.occupation_a')
ocq_b <- dbGetQuery(con, 'select SEQN, occupation from nhanes.occupation_b')
ocq_c <- dbGetQuery(con, 'select SEQN, occupation from nhanes.occupation_c')
ocq_d <- dbGetQuery(con, 'select SEQN, occupation from nhanes.occupation_d')

bt_a <- merge(bt_a, ocq_a, all.x=T)
bt_b <- merge(bt_b, ocq_b, all.x=T)
bt_c <- merge(bt_c, ocq_c, all.x=T)
bt_d <- merge(bt_d, ocq_d, all.x=T)


## recode occupation
## recode education
## recode current_past_smoking
categorical_unfold <- function(dat, varname) {
	form <- as.formula(sprintf('~factor(%s)-1', varname))
	A <- model.matrix(form,dat)
	nam <- unique(sort(as.character(dat[, varname])))
	colnames(A) <- paste(varname, nam, sep="_")
	return(A)
}

categorical_add <- function(dat, varname) {
	## adds a categorical variable as multiple binary variables to a dataset
	notNa <- !is.na(dat[, varname])
	dat.edu <- dat[notNa, c('SEQN', varname)]
	dat.edu <- cbind(dat.edu, categorical_unfold(dat.edu, varname))
	dat.edu[, varname] <- NULL
 	dat <- merge(dat, dat.edu, by.x='SEQN', by.y='SEQN', all.x=T)	
	return(dat)
}
potential_si_names <- function(colNmes) {
	lbxurx <- colNmes[grep('^LBX|^URX', colNmes)]
	siNames <- c()
	for(ii in 1:length(lbxurx)) {
		first <- substr(lbxurx[ii], 1, 2)
		ending <- substr(lbxurx[ii], 4, nchar(lbxurx[ii]))
		if(first == 'LB') {
			siNames <- c(siNames, sprintf('%sD%sSI', first, ending))
		} else {
			siNames <- c(siNames, sprintf('%sX%sSI', first, ending))
		}
		
	}
	return(siNames)
}

potential_la_names <- function(colNmes) {
	lbxurx <- colNmes[grep('^LBX', colNmes)]
	siNames <- c()
	for(ii in 1:length(lbxurx)) {
		siNames <- c(siNames, sprintf('%sLA', lbxurx[ii]))
	}
	return(siNames)
}

### filter out the variables dont need
### keep weight, SDDSRVYR, SDMVPSU, SDDSRVYR, SEQN

bt_a <- categorical_add(bt_a, 'education')
bt_b <- categorical_add(bt_b, 'education')
bt_c <- categorical_add(bt_c, 'education')
bt_d <- categorical_add(bt_d, 'education')

bt_a <- categorical_add(bt_a, 'occupation')
bt_b <- categorical_add(bt_b, 'occupation')
bt_c <- categorical_add(bt_c, 'occupation')
bt_d <- categorical_add(bt_d, 'occupation')

bt_a <- categorical_add(bt_a, 'current_past_smoking')
bt_b <- categorical_add(bt_b, 'current_past_smoking')
bt_c <- categorical_add(bt_c, 'current_past_smoking')
bt_d <- categorical_add(bt_d, 'current_past_smoking')


metaTable <- metaTable[metaTable$var != 'education', ]
metaTable <- metaTable[metaTable$var != 'current_past_smoking', ]

## now add it.
metaTable$tab_name <- NULL
metaTable$tab_desc <- NULL
metaTable$tab_desc_ewas <- NULL
metaTable$is_weight <- NULL
metaTable$version_date <- NULL
metaTable$comment_var <- NULL
metaTable$analyzable <- NULL


newRow <- function(varname, module, var_desc, var_desc_ewas,series, categorical_ref_group=NA, categorical_levels=NA, var_desc_ewas_sub=NA, binary_ref_group=NA, is_comment=0, is_questionnaire=1, is_ecological=0, is_binary=0, is_ordinal=0) {	
	data.frame(var=varname, module=module,var_desc=var_desc, series=series, var_desc_ewas=var_desc_ewas,categorical_ref_group=categorical_ref_group,categorical_levels=categorical_levels, var_desc_ewas_sub=var_desc_ewas_sub, binary_ref_group=binary_ref_group, is_comment=is_comment, is_questionnaire=is_questionnaire, is_ecological=is_ecological, is_binary=is_binary, is_ordinal=is_ordinal)
}

for(ser in c('1999-2000', '2001-2002', '2003-2004', '2005-2006')) {
	metaTable <- rbind(metaTable, newRow('education_0', 'demographics', '<HS education', 'demographics', series = ser, var_desc_ewas_sub='education'))
	metaTable <- rbind(metaTable, newRow('education_1', 'demographics', '=HS education', 'demographics',series= ser, var_desc_ewas_sub='education'))
	metaTable <- rbind(metaTable, newRow('education_2', 'demographics', '>HS education', 'demographics', series = ser, var_desc_ewas_sub='education'))

	metaTable <- rbind(metaTable, newRow('occupation_1', 'demographics', 'never worked', 'demographics',series=ser, var_desc_ewas_sub='occupation'))
	metaTable <- rbind(metaTable, newRow('occupation_2', 'demographics', 'blue-collar semi', 'demographics',series=ser,var_desc_ewas_sub='occupation'))
	metaTable <- rbind(metaTable, newRow('occupation_3', 'demographics', 'blue-collar high', 'demographics',series=ser,var_desc_ewas_sub='occupation'))
	metaTable <- rbind(metaTable, newRow('occupation_4', 'demographics', 'white-collar semi', 'demographics',series=ser, var_desc_ewas_sub='occupation'))
	metaTable <- rbind(metaTable, newRow('occupation_5', 'demographics', 'white-collar high', 'demographics',series=ser, var_desc_ewas_sub='occupation'))

	metaTable <- rbind(metaTable, newRow('current_past_smoking_0', 'questionnaire', 'never smoked', 'smoking behavior', series=ser))
	metaTable <- rbind(metaTable, newRow('current_past_smoking_1', 'questionnaire', 'past smoker', 'smoking behavior', series=ser))
	metaTable <- rbind(metaTable, newRow('current_past_smoking_2', 'questionnaire', 'current smoker', 'smoking behavior',series=ser))
}
 
allvars <- unique(c(metaTable$var, 'SEQN', 'SDMVSTRA', 'SDMVPSU', 'WTMEC2YR', 'RIDAGEYR', 'RIAGENDR', 'DMDMARTL'))
## remove the SI ones
allvars <- setdiff(allvars, potential_si_names(allvars))
allvars <- setdiff(allvars, potential_la_names(allvars))
allvars <- allvars[-grep('Unknown$', allvars)] # remove "_Unknown"
allvars <- allvars[-grep('^VTQ', allvars)] # remove "VTQ"
allvars <- allvars[-grep('^D.+Q$', allvars)] ## remove food quantity for now

bt_a <- bt_a[, colnames(bt_a) %in% allvars]
bt_b <- bt_b[, colnames(bt_b) %in% allvars]
bt_c <- bt_c[, colnames(bt_c) %in% allvars]
bt_d <- bt_d[, colnames(bt_d) %in% allvars]

### now save the data
save(bt_a, bt_b, bt_c, bt_d, allvars, metaTable, file='~/EWAS/LD/corrdata6/bigtable_correlation_nhanes.Rdata')