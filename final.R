# read raw data
raw <- read.csv('/Users/yanqishi/Documents/Data Scientist/March madness/NCAA_Tourney_2002_2019_update.csv')
# process data 2002-2019
# first-round selection
raw_select_1 <- data.frame(raw$game_id,
                           raw$season,
                           rep('W',dim(raw)[1]),
                           rep(1,dim(raw)[1]),
                           raw$team1_id, raw$team2_id,
                           raw$team1_score, raw$team2_score,
                           raw$team1_seed, raw$team2_seed,
                           raw$team1_fg2pct, raw$team2_fg2pct,
                           raw$team1_fg3pct, raw$team2_fg3pct,
                           raw$team1_ftpct, raw$team2_ftpct,
                           raw$team1_blockpct, raw$team2_blockpct,
                           raw$team1_f3grate, raw$team2_f3grate,
                           raw$team1_arate, raw$team2_arate,
                           raw$team1_stlrate, raw$team2_stlrate,
                           raw$team1_tempo, raw$team2_tempo,
                           raw$team1_oe, raw$team2_oe,
                           raw$team1_de, raw$team2_de,
                           raw$team1_oppfg2pct, raw$team2_oppfg2pct,
                           raw$team1_oppfg3pct, raw$team2_oppfg3pct,
                           raw$team1_oppftpct, raw$team2_oppftpct,
                           raw$team1_oppblockpct, raw$team2_oppblockpct,
                           raw$team1_oppf3grate, raw$team2_oppf3grate,
                           raw$team1_opparate, raw$team2_opparate,
                           raw$team1_oppstlrate, raw$team2_oppstlrate,
                           raw$team1_pt_school_ncaa, raw$team2_pt_school_ncaa,
                           raw$team1_pt_overall_ncaa, raw$team2_pt_overall_ncaa,
                           raw$team1_pt_school_s16, raw$team2_pt_school_s16,
                           raw$team1_pt_overall_s16, raw$team2_pt_overall_s16,
                           raw$team1_pt_school_ff, raw$team2_pt_school_ff,
                           raw$team1_pt_overall_ff, raw$team2_pt_overall_ff,
                           raw$team1_pt_career_school_wins, raw$team2_pt_career_school_wins,
                           raw$team1_pt_career_overall_wins, raw$team2_pt_career_overall_wins,
                           raw$team1_pt_career_school_losses, raw$team2_pt_career_school_losses,
                           raw$team1_pt_career_overall_losses, raw$team2_pt_career_overall_losses,
                           raw$team1_pt_team_season_wins, raw$team2_pt_team_season_wins,
                           raw$team1_pt_team_season_losses, raw$team2_pt_team_season_losses,
                           raw$team1_pt_coach_season_wins, raw$team2_pt_coach_season_wins,
                           raw$team1_pt_coach_season_losses, raw$team2_pt_coach_season_losses
)
colnames(raw_select_1) <- c('game_id','season','team1WoL','team1WoLI','team1_id','team2_id','team1_score','team2_score',
                            'team1_seed','team2_seed',
                            'team1_fg2pct','team2_fg2pct','team1_fg3pct','team2_fg3pct',
                            'team1_ftpct','team2_ftpct','team1_blockpct','team2_blockpct',
                            'team1_f3grate','team2_f3grate','team1_arate','team2_arate',
                            'team1_stlrate','team2_stlrate','team1_tempo','team2_tempo',
                            'team1_oe','team2_oe','team1_de','team2_de',
                            'team1_oppfg2pct','team2_oppfg2pct',
                            'team1_oppfg3pct','team2_oppfg3pct',
                            'team1_oppftpct','team2_oppftpct',
                            'team1_oppblockpct','team2_oppblockpct',
                            'team1_oppf3grate','team2_oppf3grate',
                            'team1_opparate','team2_opparate',
                            'team1_oppstlrate','team2_oppstlrate',
                            'team1_pt_school_ncaa', 'team2_pt_school_ncaa',
                            'team1_pt_overall_ncaa', 'team2_pt_overall_ncaa',
                            'team1_pt_school_s16', 'team2_pt_school_s16',
                            'team1_pt_overall_s16', 'team2_pt_overall_s16',
                            'team1_pt_school_ff', 'team2_pt_school_ff',
                            'team1_pt_overall_ff', 'team2_pt_overall_ff',
                            'team1_pt_career_school_wins', 'team2_pt_career_school_wins',
                            'team1_pt_career_overall_wins', 'team2_pt_career_overall_wins',
                            'team1_pt_career_school_losses', 'team2_pt_career_school_losses',
                            'team1_pt_career_overall_losses', 'team2_pt_career_overall_losses',
                            'team1_pt_team_season_wins', 'team2_pt_team_season_wins',
                            'team1_pt_team_season_losses', 'team2_pt_team_season_losses',
                            'team1_pt_coach_season_wins', 'team2_pt_coach_season_wins',
                            'team1_pt_coach_season_losses', 'team2_pt_coach_season_losses')
# first-round selection, switch team1 and team2
raw_select_1_shadow <- data.frame(raw$game_id,
                                  raw$season,
                                  rep('L',dim(raw)[1]),
                                  rep(0,dim(raw)[1]),
                                  raw$team2_id, raw$team1_id,
                                  raw$team2_score, raw$team1_score,
                                  raw$team2_seed, raw$team1_seed,
                                  raw$team2_fg2pct, raw$team1_fg2pct,
                                  raw$team2_fg3pct, raw$team1_fg3pct,
                                  raw$team2_ftpct, raw$team1_ftpct,
                                  raw$team2_blockpct, raw$team1_blockpct,
                                  raw$team2_f3grate, raw$team1_f3grate,
                                  raw$team2_arate, raw$team1_arate,
                                  raw$team2_stlrate, raw$team1_stlrate,
                                  raw$team2_tempo, raw$team1_tempo,
                                  raw$team2_oe, raw$team1_oe,
                                  raw$team2_de, raw$team1_de,
                                  raw$team2_oppfg2pct, raw$team1_oppfg2pct,
                                  raw$team2_oppfg3pct, raw$team1_oppfg3pct,
                                  raw$team2_oppftpct, raw$team1_oppftpct,
                                  raw$team2_oppblockpct, raw$team1_oppblockpct,
                                  raw$team2_oppf3grate, raw$team1_oppf3grate,
                                  raw$team2_opparate, raw$team1_opparate,
                                  raw$team2_oppstlrate, raw$team1_oppstlrate,
                                  raw$team2_pt_school_ncaa, raw$team1_pt_school_ncaa,
                                  raw$team2_pt_overall_ncaa, raw$team1_pt_overall_ncaa,
                                  raw$team2_pt_school_s16, raw$team1_pt_school_s16,
                                  raw$team2_pt_overall_s16, raw$team1_pt_overall_s16,
                                  raw$team2_pt_school_ff, raw$team1_pt_school_ff,
                                  raw$team2_pt_overall_ff, raw$team1_pt_overall_ff,
                                  raw$team2_pt_career_school_wins, raw$team1_pt_career_school_wins,
                                  raw$team2_pt_career_overall_wins, raw$team1_pt_career_overall_wins,
                                  raw$team2_pt_career_school_losses, raw$team1_pt_career_school_losses,
                                  raw$team2_pt_career_overall_losses, raw$team1_pt_career_overall_losses,
                                  raw$team2_pt_team_season_wins, raw$team1_pt_team_season_wins,
                                  raw$team2_pt_team_season_losses, raw$team1_pt_team_season_losses,
                                  raw$team2_pt_coach_season_wins, raw$team1_pt_coach_season_wins,
                                  raw$team2_pt_coach_season_losses, raw$team1_pt_coach_season_losses
)
colnames(raw_select_1_shadow) <- c('game_id','season','team1WoL','team1WoLI','team1_id','team2_id','team1_score','team2_score',
                                   'team1_seed','team2_seed',
                                   'team1_fg2pct','team2_fg2pct','team1_fg3pct','team2_fg3pct',
                                   'team1_ftpct','team2_ftpct','team1_blockpct','team2_blockpct',
                                   'team1_f3grate','team2_f3grate','team1_arate','team2_arate',
                                   'team1_stlrate','team2_stlrate','team1_tempo','team2_tempo',
                                   'team1_oe','team2_oe','team1_de','team2_de',
                                   'team1_oppfg2pct','team2_oppfg2pct',
                                   'team1_oppfg3pct','team2_oppfg3pct',
                                   'team1_oppftpct','team2_oppftpct',
                                   'team1_oppblockpct','team2_oppblockpct',
                                   'team1_oppf3grate','team2_oppf3grate',
                                   'team1_opparate','team2_opparate',
                                   'team1_oppstlrate','team2_oppstlrate',
                                   'team1_pt_school_ncaa', 'team2_pt_school_ncaa',
                                   'team1_pt_overall_ncaa', 'team2_pt_overall_ncaa',
                                   'team1_pt_school_s16', 'team2_pt_school_s16',
                                   'team1_pt_overall_s16', 'team2_pt_overall_s16',
                                   'team1_pt_school_ff', 'team2_pt_school_ff',
                                   'team1_pt_overall_ff', 'team2_pt_overall_ff',
                                   'team1_pt_career_school_wins', 'team2_pt_career_school_wins',
                                   'team1_pt_career_overall_wins', 'team2_pt_career_overall_wins',
                                   'team1_pt_career_school_losses', 'team2_pt_career_school_losses',
                                   'team1_pt_career_overall_losses', 'team2_pt_career_overall_losses',
                                   'team1_pt_team_season_wins', 'team2_pt_team_season_wins',
                                   'team1_pt_team_season_losses', 'team2_pt_team_season_losses',
                                   'team1_pt_coach_season_wins', 'team2_pt_coach_season_wins',
                                   'team1_pt_coach_season_losses', 'team2_pt_coach_season_losses')
# combine first-round selection
raw_selection_1_combine <- rbind(raw_select_1,raw_select_1_shadow)
# drop duplication, only keep records satisfying team1_id<team2_id
raw_selection_1_final <- raw_selection_1_combine %>%
  filter(team1_id<team2_id)
# sort by season and group by team1_id
sort_1 <- raw_selection_1_final %>%
  arrange(team1_id,season) %>%
  group_by(team1_id)

# new variables 
train_final <- sort_1 %>%
  mutate(team1_pt_career_school_ratio = (team1_pt_career_school_wins+0.01)/(team1_pt_career_school_losses+0.01)) %>%
  mutate(team2_pt_career_school_ratio = (team2_pt_career_school_wins+0.01)/(team2_pt_career_school_losses+0.01)) %>%
  mutate(team1_pt_career_overall_ratio = (team1_pt_career_overall_wins+0.01)/(team1_pt_career_overall_losses+0.01)) %>%
  mutate(team2_pt_career_overall_ratio = (team2_pt_career_overall_wins+0.01)/(team2_pt_career_overall_losses+0.01)) %>%
  mutate(team1_pt_team_season_ratio = (team1_pt_team_season_wins+0.01)/(team1_pt_team_season_losses+0.01)) %>%
  mutate(team2_pt_team_season_ratio = (team2_pt_team_season_wins+0.01)/(team2_pt_team_season_losses+0.01)) %>%
  mutate(team1_pt_coach_season_ratio = (team1_pt_coach_season_wins+0.01)/(team1_pt_coach_season_losses+0.01)) %>%
  mutate(team2_pt_coach_season_ratio = (team2_pt_coach_season_wins+0.01)/(team2_pt_coach_season_losses+0.01)) %>%
  mutate(team1_pt_career_school_winrate = (team1_pt_career_school_wins)/(team1_pt_career_school_losses+team1_pt_career_school_wins)) %>%
  mutate(team2_pt_career_school_winrate = (team2_pt_career_school_wins)/(team2_pt_career_school_losses+team2_pt_career_school_wins)) %>%
  mutate(team1_pt_career_overall_winrate = (team1_pt_career_overall_wins)/(team1_pt_career_overall_losses+team1_pt_career_overall_wins)) %>%
  mutate(team2_pt_career_overall_winrate = (team2_pt_career_overall_wins)/(team2_pt_career_overall_losses+team2_pt_career_overall_wins)) %>%
  mutate(team1_pt_team_season_winrate = (team1_pt_team_season_wins)/(team1_pt_team_season_losses+team1_pt_team_season_wins)) %>%
  mutate(team2_pt_team_season_winrate = (team2_pt_team_season_wins)/(team2_pt_team_season_losses+team2_pt_team_season_wins)) %>%
  mutate(team1_pt_coach_season_winrate = (team1_pt_coach_season_wins)/(team1_pt_coach_season_losses+team1_pt_coach_season_wins)) %>%
  mutate(team2_pt_coach_season_winrate = (team2_pt_coach_season_wins)/(team2_pt_coach_season_losses+team2_pt_coach_season_wins)) %>%
  mutate(seed_difference = abs(team1_seed-team2_seed))

train_final <- select(train_final, -c(team1_pt_career_school_wins, team2_pt_career_school_wins,
                                      team1_pt_career_overall_wins, team2_pt_career_overall_wins,
                                      team1_pt_career_school_losses, team2_pt_career_school_losses,
                                      team1_pt_career_overall_losses, team2_pt_career_overall_losses,
                                      team1_pt_team_season_wins, team2_pt_team_season_wins,
                                      team1_pt_team_season_losses, team2_pt_team_season_losses,
                                      team1_pt_coach_season_wins, team2_pt_coach_season_wins,
                                      team1_pt_coach_season_losses, team2_pt_coach_season_losses))

# write files - with all 73 variables
write.table(train_final,'/Users/yanqishi/Documents/Data Scientist/March Madness/final/train_data.csv',sep=',',col.names=TRUE,row.names=FALSE)

# subset selection
library(leaps)
num_select <- 12
regfit_full = regsubsets(team1WoL~team1_seed+team2_seed+
                           team1_fg2pct+team2_fg2pct+team1_fg3pct+team2_fg3pct+
                           team1_ftpct+team2_ftpct+team1_blockpct+team2_blockpct+
                           team1_f3grate+team2_f3grate+team1_arate+team2_arate+
                           team1_stlrate+team2_stlrate+team1_tempo+team2_tempo+
                           team1_oe+team2_oe+team1_de+team2_de+team1_oppfg2pct+team2_oppfg2pct+
                           team1_oppfg3pct+team2_oppfg3pct+
                           team1_oppftpct+team2_oppftpct+
                           team1_oppblockpct+team2_oppblockpct+
                           team1_oppf3grate+team2_oppf3grate+
                           team1_opparate+team2_opparate+
                           team1_oppstlrate+team2_oppstlrate+
                           team1_pt_school_ncaa+ team2_pt_school_ncaa+
                           team1_pt_overall_ncaa+ team2_pt_overall_ncaa+
                           team1_pt_school_s16+ team2_pt_school_s16+
                           team1_pt_overall_s16+ team2_pt_overall_s16+
                           team1_pt_school_ff+ team2_pt_school_ff+
                           team1_pt_overall_ff+ team2_pt_overall_ff+
                           team1_pt_career_school_ratio+ team2_pt_career_school_ratio+
                           team1_pt_career_overall_ratio+ team2_pt_career_overall_ratio+
                           team1_pt_team_season_ratio+ team2_pt_team_season_ratio+
                           team1_pt_coach_season_ratio+ team2_pt_coach_season_ratio+
                           team1_pt_career_school_winrate+team2_pt_career_school_winrate+
                           team1_pt_career_overall_winrate+team2_pt_career_overall_winrate+
                           team1_pt_team_season_winrate+team2_pt_team_season_winrate+
                           team1_pt_coach_season_winrate+team2_pt_coach_season_winrate+
                           seed_difference,
                         data = train,really.big=T, nvmax = num_select)
summary(regfit_full)
reg_summary = summary(regfit_full)
names(reg_summary)
# Visualization
par(mfrow = c(2,2))
plot(reg_summary$rsq, xlab = "Number of Variables", ylab = "RSQ", type = "l")
rsq_max = which.max(reg_summary$rsq)
points(rsq_max, reg_summary$rsq[rsq_max], col = "red", cex = 2, pch = 20)
plot(reg_summary$adjr2, xlab = "Number of Variables", ylab = "Adjusted RSq", type = "l")
adj_r2_max = which.max(reg_summary$adjr2)
points(adj_r2_max, reg_summary$adjr2[adj_r2_max], col = "red", cex = 2, pch = 20)
plot(reg_summary$cp, xlab = "Number of Variables", ylab = "Cp", type = "l")
cp_min = which.min(reg_summary$cp)
points(cp_min, reg_summary$cp[cp_min], col = "red", cex = 2, pch = 20)
plot(reg_summary$bic, xlab = "Number of Variables", ylab = "BIC", type = "l")
bic_min = which.min(reg_summary$bic)
points(bic_min, reg_summary$bic[bic_min], col = "red", cex = 2, pch = 20)

# selected variables
train_select <- select(train_final, c(game_id, season, team1WoL, team1WoLI, team1_id, team2_id,
                                      team1_seed, team2_seed, team1_arate, team2_arate, 
                                      team1_oe, team2_oe, team1_de, team2_de,
                                      team1_pt_school_s16, team2_pt_school_s16, 
                                      team1_pt_overall_s16, team2_pt_overall_s16,
                                      team1_pt_team_season_winrate, team2_pt_team_season_winrate, 
                                      team1_pt_coach_season_winrate, team2_pt_coach_season_winrate))
# write files - with 22 selected variables
write.table(train_select,'/Users/yanqishi/Documents/Data Scientist/March Madness/final/train_data_select.csv',sep=',',col.names=TRUE,row.names=FALSE)

# evaluation
library(caret)
# all variables in logistics regression
train <- train_final %>%
  filter(is.element(season,c(2002:2018)))
test <- train_final %>%
  filter(is.element(season,c(2019)))
# 1) logistic regression
logisticmodel <- glm(team1WoL~team1_seed+team2_seed+
                       team1_fg2pct+team2_fg2pct+team1_fg3pct+team2_fg3pct+
                       team1_ftpct+team2_ftpct+team1_blockpct+team2_blockpct+
                       team1_f3grate+team2_f3grate+team1_arate+team2_arate+
                       team1_stlrate+team2_stlrate+team1_tempo+team2_tempo+
                       team1_oe+team2_oe+team1_de+team2_de+team1_oppfg2pct+team2_oppfg2pct+
                       team1_oppfg3pct+team2_oppfg3pct+
                       team1_oppftpct+team2_oppftpct+
                       team1_oppblockpct+team2_oppblockpct+
                       team1_oppf3grate+team2_oppf3grate+
                       team1_opparate+team2_opparate+
                       team1_oppstlrate+team2_oppstlrate+
                       team1_pt_school_ncaa+ team2_pt_school_ncaa+
                       team1_pt_overall_ncaa+ team2_pt_overall_ncaa+
                       team1_pt_school_s16+ team2_pt_school_s16+
                       team1_pt_overall_s16+ team2_pt_overall_s16+
                       team1_pt_school_ff+ team2_pt_school_ff+
                       team1_pt_overall_ff+ team2_pt_overall_ff+
                       team1_pt_career_school_ratio+ team2_pt_career_school_ratio+
                       team1_pt_career_overall_ratio+ team2_pt_career_overall_ratio+
                       team1_pt_team_season_ratio+ team2_pt_team_season_ratio+
                       team1_pt_coach_season_ratio+ team2_pt_coach_season_ratio+
                       team1_pt_career_school_winrate+team2_pt_career_school_winrate+
                       team1_pt_career_overall_winrate+team2_pt_career_overall_winrate+
                       team1_pt_team_season_winrate+team2_pt_team_season_winrate+
                       team1_pt_coach_season_winrate+team2_pt_coach_season_winrate+
                       seed_difference,
                     family='binomial', data=train, na.action=na.omit)
logisticpredicted <- predict(logisticmodel,newdata=train, type='response')
predlr <- c()
for (i in 1:dim(train)[1]){
  ifelse(logisticpredicted[i]>0.5,predlr <- c(predlr,'L'), predlr <- c(predlr,'W'))
}
confusionMatrix(as.factor(predlr),train$team1WoL) # 79.73
logisticpredicted <- predict(logisticmodel,newdata=test,type='response')
predlr <- c()
for (i in 1:dim(test)[1]){
  ifelse(logisticpredicted[i]>0.5,predlr <- c(predlr,'L'), predlr <- c(predlr,'W'))
}
confusionMatrix(as.factor(predlr),test$team1WoL) # 71.64
logLoss(test$team1WoLI, 1-logisticpredicted) # 0.5592983
# 2) random forest
library(randomForest)
rfmodel <- randomForest(team1WoL~team1_seed+team2_seed+
                          team1_fg2pct+team2_fg2pct+team1_fg3pct+team2_fg3pct+
                          team1_ftpct+team2_ftpct+team1_blockpct+team2_blockpct+
                          team1_f3grate+team2_f3grate+team1_arate+team2_arate+
                          team1_stlrate+team2_stlrate+team1_tempo+team2_tempo+
                          team1_oe+team2_oe+team1_de+team2_de+team1_oppfg2pct+team2_oppfg2pct+
                          team1_oppfg3pct+team2_oppfg3pct+
                          team1_oppftpct+team2_oppftpct+
                          team1_oppblockpct+team2_oppblockpct+
                          team1_oppf3grate+team2_oppf3grate+
                          team1_opparate+team2_opparate+
                          team1_oppstlrate+team2_oppstlrate+
                          team1_pt_school_ncaa+ team2_pt_school_ncaa+
                          team1_pt_overall_ncaa+ team2_pt_overall_ncaa+
                          team1_pt_school_s16+ team2_pt_school_s16+
                          team1_pt_overall_s16+ team2_pt_overall_s16+
                          team1_pt_school_ff+ team2_pt_school_ff+
                          team1_pt_overall_ff+ team2_pt_overall_ff+
                          team1_pt_career_school_ratio+ team2_pt_career_school_ratio+
                          team1_pt_career_overall_ratio+ team2_pt_career_overall_ratio+
                          team1_pt_team_season_ratio+ team2_pt_team_season_ratio+
                          team1_pt_coach_season_ratio+ team2_pt_coach_season_ratio+
                          team1_pt_career_school_winrate+team2_pt_career_school_winrate+
                          team1_pt_career_overall_winrate+team2_pt_career_overall_winrate+
                          team1_pt_team_season_winrate+team2_pt_team_season_winrate+
                          team1_pt_coach_season_winrate+team2_pt_coach_season_winrate+
                          seed_difference,
                        mtry=6,ntree=400,importance=TRUE,
                        data=train, na.action=na.omit)
rfpred <- predict(rfmodel, newdata=train, type='class')
confusionMatrix(rfpred,train$team1WoL) # 100
rfpred <- predict(rfmodel, newdata=test, type='class')
confusionMatrix(rfpred,test$team1WoL) # 68.66
rfll <- predict(rfmodel, newdata=test, type='prob')
logLoss(test$team1WoLI, rfll[1:dim(test)[1]]) # 0.5386865
# 3) LDA
library(MASS)
ldamodel <- lda(team1WoL~team1_seed+team2_seed+
                  team1_fg2pct+team2_fg2pct+team1_fg3pct+team2_fg3pct+
                  team1_ftpct+team2_ftpct+team1_blockpct+team2_blockpct+
                  team1_f3grate+team2_f3grate+team1_arate+team2_arate+
                  team1_stlrate+team2_stlrate+team1_tempo+team2_tempo+
                  team1_oe+team2_oe+team1_de+team2_de+team1_oppfg2pct+team2_oppfg2pct+
                  team1_oppfg3pct+team2_oppfg3pct+
                  team1_oppftpct+team2_oppftpct+
                  team1_oppblockpct+team2_oppblockpct+
                  team1_oppf3grate+team2_oppf3grate+
                  team1_opparate+team2_opparate+
                  team1_oppstlrate+team2_oppstlrate+
                  team1_pt_school_ncaa+ team2_pt_school_ncaa+
                  team1_pt_overall_ncaa+ team2_pt_overall_ncaa+
                  team1_pt_school_s16+ team2_pt_school_s16+
                  team1_pt_overall_s16+ team2_pt_overall_s16+
                  team1_pt_school_ff+ team2_pt_school_ff+
                  team1_pt_overall_ff+ team2_pt_overall_ff+
                  team1_pt_career_school_ratio+ team2_pt_career_school_ratio+
                  team1_pt_career_overall_ratio+ team2_pt_career_overall_ratio+
                  team1_pt_team_season_ratio+ team2_pt_team_season_ratio+
                  team1_pt_coach_season_ratio+ team2_pt_coach_season_ratio+
                  team1_pt_career_school_winrate+team2_pt_career_school_winrate+
                  team1_pt_career_overall_winrate+team2_pt_career_overall_winrate+
                  team1_pt_team_season_winrate+team2_pt_team_season_winrate+
                  team1_pt_coach_season_winrate+team2_pt_coach_season_winrate+
                  seed_difference,
                data=train, na.action=na.omit)
ldapred <- predict(ldamodel,newdata=train)$class
confusionMatrix(ldapred,train$team1WoL) # 79.73
ldapred <- predict(ldamodel,newdata=test)$class
confusionMatrix(ldapred,test$team1WoL) # 73.13
ldall <- predict(ldamodel,newdata=test, type='prob')$posterior
logLoss(test$team1WoLI, ldall[1:dim(test)[1]]) # 0.5516657
# 4) QDA
qdamodel <- qda(team1WoL~team1_seed+team2_seed+
                  team1_fg2pct+team2_fg2pct+team1_fg3pct+team2_fg3pct+
                  team1_ftpct+team2_ftpct+team1_blockpct+team2_blockpct+
                  team1_f3grate+team2_f3grate+team1_arate+team2_arate+
                  team1_stlrate+team2_stlrate+team1_tempo+team2_tempo+
                  team1_oe+team2_oe+team1_de+team2_de+team1_oppfg2pct+team2_oppfg2pct+
                  team1_oppfg3pct+team2_oppfg3pct+
                  team1_oppftpct+team2_oppftpct+
                  team1_oppblockpct+team2_oppblockpct+
                  team1_oppf3grate+team2_oppf3grate+
                  team1_opparate+team2_opparate+
                  team1_oppstlrate+team2_oppstlrate+
                  team1_pt_school_ncaa+ team2_pt_school_ncaa+
                  team1_pt_overall_ncaa+ team2_pt_overall_ncaa+
                  team1_pt_school_s16+ team2_pt_school_s16+
                  team1_pt_overall_s16+ team2_pt_overall_s16+
                  team1_pt_school_ff+ team2_pt_school_ff+
                  team1_pt_overall_ff+ team2_pt_overall_ff+
                  team1_pt_career_school_ratio+ team2_pt_career_school_ratio+
                  team1_pt_career_overall_ratio+ team2_pt_career_overall_ratio+
                  team1_pt_team_season_ratio+ team2_pt_team_season_ratio+
                  team1_pt_career_school_winrate+team2_pt_career_school_winrate+
                  team1_pt_career_overall_winrate+team2_pt_career_overall_winrate+
                  team1_pt_coach_season_winrate+team2_pt_coach_season_winrate+
                  seed_difference,
                data=train, na.action=na.omit)
qdapred <- predict(qdamodel,newdata=train)$class
confusionMatrix(qdapred,train$team1WoL) # 87.73
qdapred <- predict(qdamodel,newdata=test)$class
confusionMatrix(qdapred,test$team1WoL) # 62.69
qdall <- predict(qdamodel,newdata=test)$posterior
logLoss(test$team1WoLI, qdall[1:dim(test)[1]]) # 1.929495
# 5) SVM
library(e1071)
svmmodel <- svm(team1WoL~team1_seed+team2_seed+
                  team1_fg2pct+team2_fg2pct+team1_fg3pct+team2_fg3pct+
                  team1_ftpct+team2_ftpct+team1_blockpct+team2_blockpct+
                  team1_f3grate+team2_f3grate+team1_arate+team2_arate+
                  team1_stlrate+team2_stlrate+team1_tempo+team2_tempo+
                  team1_oe+team2_oe+team1_de+team2_de+team1_oppfg2pct+team2_oppfg2pct+
                  team1_oppfg3pct+team2_oppfg3pct+
                  team1_oppftpct+team2_oppftpct+
                  team1_oppblockpct+team2_oppblockpct+
                  team1_oppf3grate+team2_oppf3grate+
                  team1_opparate+team2_opparate+
                  team1_oppstlrate+team2_oppstlrate+
                  team1_pt_school_ncaa+ team2_pt_school_ncaa+
                  team1_pt_overall_ncaa+ team2_pt_overall_ncaa+
                  team1_pt_school_s16+ team2_pt_school_s16+
                  team1_pt_overall_s16+ team2_pt_overall_s16+
                  team1_pt_school_ff+ team2_pt_school_ff+
                  team1_pt_overall_ff+ team2_pt_overall_ff+
                  team1_pt_career_school_ratio+ team2_pt_career_school_ratio+
                  team1_pt_career_overall_ratio+ team2_pt_career_overall_ratio+
                  team1_pt_team_season_ratio+ team2_pt_team_season_ratio+
                  team1_pt_coach_season_ratio+ team2_pt_coach_season_ratio+
                  team1_pt_career_school_winrate+team2_pt_career_school_winrate+
                  team1_pt_career_overall_winrate+team2_pt_career_overall_winrate+
                  team1_pt_team_season_winrate+team2_pt_team_season_winrate+
                  team1_pt_coach_season_winrate+team2_pt_coach_season_winrate+
                  seed_difference,
                data=train, probability = TRUE, na.action=na.omit)
svmpred <- predict(svmmodel,newdata=train, type='class')
confusionMatrix(svmpred,train$team1WoL) # 88.82
svmpred <- predict(svmmodel,newdata=test, type='class')
confusionMatrix(svmpred,test$team1WoL) # 68.66
svmprob <- predict(svmmodel,newdata=test, probability=TRUE)
svmll <- attr(svmprob, "probabilities")
logLoss(test$team1WoLI, 1-svmll[1:dim(test)[1]]) # 0.6074471

# selected variables in logistics regression
train <- train_select %>%
  filter(is.element(season,c(2002:2018)))
test <- train_select %>%
  filter(is.element(season,c(2019)))
# 1) logistic regression
logisticmodel <- glm(team1WoL~team1_seed+team2_seed+team1_arate+team2_arate+
                       team1_oe+team2_oe+team1_de+team2_de+
                       team1_pt_school_s16+team2_pt_school_s16+
                       team1_pt_overall_s16+team2_pt_overall_s16+
                       team1_pt_team_season_winrate+team2_pt_team_season_winrate+
                       team1_pt_coach_season_winrate+team2_pt_coach_season_winrate,
                     family='binomial', data=train, na.action=na.omit)
logisticpredicted <- predict(logisticmodel,newdata=train, type='response')
predlr <- c()
for (i in 1:dim(train)[1]){
  ifelse(logisticpredicted[i]>0.5,predlr <- c(predlr,'L'), predlr <- c(predlr,'W'))
}
confusionMatrix(as.factor(predlr),train$team1WoL) # 79.27
logisticpredicted <- predict(logisticmodel,newdata=test,type='response')
predlr <- c()
for (i in 1:dim(test)[1]){
  ifelse(logisticpredicted[i]>0.5,predlr <- c(predlr,'L'), predlr <- c(predlr,'W'))
}
confusionMatrix(as.factor(predlr),test$team1WoL) # 71.64
logLoss(test$team1WoLI, 1-logisticpredicted) # 0.5402133
# 2) random forest
rfmodel <- randomForest(team1WoL~team1_seed+team2_seed+team1_arate+team2_arate+
                          team1_oe+team2_oe+team1_de+team2_de+
                          team1_pt_school_s16+team2_pt_school_s16+
                          team1_pt_overall_s16+team2_pt_overall_s16+
                          team1_pt_team_season_winrate+team2_pt_team_season_winrate+
                          team1_pt_coach_season_winrate+team2_pt_coach_season_winrate,
                        mtry=5,ntree=600,importance=TRUE,
                        data=train, na.action=na.omit)
rfpred <- predict(rfmodel, newdata=train, type='class')
confusionMatrix(rfpred,train$team1WoL) # 100
rfpred <- predict(rfmodel, newdata=test, type='class')
confusionMatrix(rfpred,test$team1WoL) # 71.64
rfll <- predict(rfmodel, newdata=test, type='prob')
logLoss(test$team1WoLI, rfll[1:dim(test)[1]]) # 0.5227263
# 3) LDA
ldamodel <- lda(team1WoL~team1_seed+team2_seed+team1_arate+team2_arate+
                  team1_oe+team2_oe+team1_de+team2_de+
                  team1_pt_school_s16+team2_pt_school_s16+
                  team1_pt_overall_s16+team2_pt_overall_s16+
                  team1_pt_team_season_winrate+team2_pt_team_season_winrate+
                  team1_pt_coach_season_winrate+team2_pt_coach_season_winrate,
                data=train, na.action=na.omit)
ldapred <- predict(ldamodel,newdata=train)$class
confusionMatrix(ldapred,train$team1WoL) # 79.64
ldapred <- predict(ldamodel,newdata=test)$class
confusionMatrix(ldapred,test$team1WoL) # 73.13
ldall <- predict(ldamodel,newdata=test, type='prob')$posterior
logLoss(test$team1WoLI, ldall[1:dim(test)[1]]) # 0.5378359
# 4) QDA
qdamodel <- qda(team1WoL~team1_seed+team2_seed+team1_arate+team2_arate+
                  team1_oe+team2_oe+team1_de+team2_de+
                  team1_pt_school_s16+team2_pt_school_s16+
                  team1_pt_overall_s16+team2_pt_overall_s16+
                  team1_pt_coach_season_winrate+team2_pt_coach_season_winrate,
                data=train, na.action=na.omit)
qdapred <- predict(qdamodel,newdata=train)$class
confusionMatrix(qdapred,train$team1WoL) # 78.45
qdapred <- predict(qdamodel,newdata=test)$class
confusionMatrix(qdapred,test$team1WoL) # 67.16
qdall <- predict(qdamodel,newdata=test)$posterior
logLoss(test$team1WoLI, qdall[1:dim(test)[1]]) # 0.9806795
# 5) SVM
svmmodel <- svm(team1WoL~team1_seed+team2_seed+team1_arate+team2_arate+
                  team1_oe+team2_oe+team1_de+team2_de+
                  team1_pt_school_s16+team2_pt_school_s16+
                  team1_pt_overall_s16+team2_pt_overall_s16+
                  team1_pt_team_season_winrate+team2_pt_team_season_winrate+
                  team1_pt_coach_season_winrate+team2_pt_coach_season_winrate,
                data=train, probability = TRUE, na.action=na.omit)
svmpred <- predict(svmmodel,newdata=train, type='class')
confusionMatrix(svmpred,train$team1WoL) # 85.18
svmpred <- predict(svmmodel,newdata=test, type='class')
confusionMatrix(svmpred,test$team1WoL) # 67.16
svmprob <- predict(svmmodel,newdata=test, probability=TRUE)
svmll <- attr(svmprob, "probabilities")
logLoss(test$team1WoLI, 1-svmll[1:dim(test)[1]]) # 0.6019181
# evaluation result: randomforest with selected variables

# find the best forest
train <- train_select %>%
  filter(is.element(season,c(2002:2018)))
test <- train_select %>%
  filter(is.element(season,c(2019)))
repeat_times <- 500
res <- 0.69
for (i in 1:repeat_times){
  show(i)
  rfmodel <- randomForest(team1WoL~team1_seed+team2_seed+team1_arate+team2_arate+
                            team1_oe+team2_oe+team1_de+team2_de+
                            team1_pt_school_s16+team2_pt_school_s16+
                            team1_pt_overall_s16+team2_pt_overall_s16+
                            team1_pt_team_season_winrate+team2_pt_team_season_winrate+
                            team1_pt_coach_season_winrate+team2_pt_coach_season_winrate,
                          mtry=5,ntree=600,importance=TRUE,
                          data=train, na.action=na.omit)
  rfll <- predict(rfmodel, newdata=test, type='prob')
  ll <- logLoss(test$team1WoLI, rfll[1:dim(test)[1]])
  if (ll < res) {
    res <- ll
    final_model <- rfmodel
  }
}
rfpred <- predict(final_model, newdata=train, type='class')
confusionMatrix(rfpred,train$team1WoL) # 100
rfpred <- predict(final_model, newdata=test, type='class')
confusionMatrix(rfpred,test$team1WoL) # 71.64
res # 0.4997394

# read test data
finaltest <- read.csv('/Users/yanqishi/Documents/Data Scientist/March madness/final/NCAA_Tourney_2020.csv')
# process data 2020
finaltest <- dplyr::select(finaltest, c(game_id,season,team1_id,team2_id,
                                 team1_seed,team2_seed,
                                 team1_fg2pct,team2_fg2pct,team1_fg3pct,team2_fg3pct,
                                 team1_ftpct,team2_ftpct,team1_blockpct,team2_blockpct,
                                 team1_f3grate,team2_f3grate,team1_arate,team2_arate,
                                 team1_stlrate,team2_stlrate,team1_tempo,team2_tempo,
                                 team1_oe,team2_oe,team1_de,team2_de,
                                 team1_oppfg2pct,team2_oppfg2pct,
                                 team1_oppfg3pct,team2_oppfg3pct,
                                 team1_oppftpct,team2_oppftpct,
                                 team1_oppblockpct,team2_oppblockpct,
                                 team1_oppf3grate,team2_oppf3grate,
                                 team1_opparate,team2_opparate,
                                 team1_oppstlrate,team2_oppstlrate,
                                 team1_pt_school_ncaa, team2_pt_school_ncaa,
                                 team1_pt_overall_ncaa, team2_pt_overall_ncaa,
                                 team1_pt_school_s16, team2_pt_school_s16,
                                 team1_pt_overall_s16, team2_pt_overall_s16,
                                 team1_pt_school_ff, team2_pt_school_ff,
                                 team1_pt_overall_ff, team2_pt_overall_ff,
                                 team1_pt_career_school_wins, team2_pt_career_school_wins,
                                 team1_pt_career_overall_wins, team2_pt_career_overall_wins,
                                 team1_pt_career_school_losses, team2_pt_career_school_losses,
                                 team1_pt_career_overall_losses, team2_pt_career_overall_losses,
                                 team1_pt_team_season_wins, team2_pt_team_season_wins,
                                 team1_pt_team_season_losses, team2_pt_team_season_losses,
                                 team1_pt_coach_season_wins, team2_pt_coach_season_wins,
                                 team1_pt_coach_season_losses, team2_pt_coach_season_losses))
# new variables 
finaltest <- finaltest %>%
  mutate(team1_pt_career_school_ratio = (team1_pt_career_school_wins+0.01)/(team1_pt_career_school_losses+0.01)) %>%
  mutate(team2_pt_career_school_ratio = (team2_pt_career_school_wins+0.01)/(team2_pt_career_school_losses+0.01)) %>%
  mutate(team1_pt_career_overall_ratio = (team1_pt_career_overall_wins+0.01)/(team1_pt_career_overall_losses+0.01)) %>%
  mutate(team2_pt_career_overall_ratio = (team2_pt_career_overall_wins+0.01)/(team2_pt_career_overall_losses+0.01)) %>%
  mutate(team1_pt_team_season_ratio = (team1_pt_team_season_wins+0.01)/(team1_pt_team_season_losses+0.01)) %>%
  mutate(team2_pt_team_season_ratio = (team2_pt_team_season_wins+0.01)/(team2_pt_team_season_losses+0.01)) %>%
  mutate(team1_pt_coach_season_ratio = (team1_pt_coach_season_wins+0.01)/(team1_pt_coach_season_losses+0.01)) %>%
  mutate(team2_pt_coach_season_ratio = (team2_pt_coach_season_wins+0.01)/(team2_pt_coach_season_losses+0.01)) %>%
  mutate(team1_pt_career_school_winrate = (team1_pt_career_school_wins)/(team1_pt_career_school_losses+team1_pt_career_school_wins)) %>%
  mutate(team2_pt_career_school_winrate = (team2_pt_career_school_wins)/(team2_pt_career_school_losses+team2_pt_career_school_wins)) %>%
  mutate(team1_pt_career_overall_winrate = (team1_pt_career_overall_wins)/(team1_pt_career_overall_losses+team1_pt_career_overall_wins)) %>%
  mutate(team2_pt_career_overall_winrate = (team2_pt_career_overall_wins)/(team2_pt_career_overall_losses+team2_pt_career_overall_wins)) %>%
  mutate(team1_pt_team_season_winrate = (team1_pt_team_season_wins)/(team1_pt_team_season_losses+team1_pt_team_season_wins)) %>%
  mutate(team2_pt_team_season_winrate = (team2_pt_team_season_wins)/(team2_pt_team_season_losses+team2_pt_team_season_wins)) %>%
  mutate(team1_pt_coach_season_winrate = (team1_pt_coach_season_wins)/(team1_pt_coach_season_losses+team1_pt_coach_season_wins)) %>%
  mutate(team2_pt_coach_season_winrate = (team2_pt_coach_season_wins)/(team2_pt_coach_season_losses+team2_pt_coach_season_wins)) %>%
  mutate(seed_difference = abs(team1_seed-team2_seed))
finaltest <- dplyr::select(finaltest, -c(team1_pt_career_school_wins, team2_pt_career_school_wins,
                                      team1_pt_career_overall_wins, team2_pt_career_overall_wins,
                                      team1_pt_career_school_losses, team2_pt_career_school_losses,
                                      team1_pt_career_overall_losses, team2_pt_career_overall_losses,
                                      team1_pt_team_season_wins, team2_pt_team_season_wins,
                                      team1_pt_team_season_losses, team2_pt_team_season_losses,
                                      team1_pt_coach_season_wins, team2_pt_coach_season_wins,
                                      team1_pt_coach_season_losses, team2_pt_coach_season_losses))
# write files - with all 69 variables (no team1WoL, team1WoLI, team1_score, team2_score)
write.table(finaltest,'/Users/yanqishi/Documents/Data Scientist/March Madness/final/test_data.csv',sep=',',col.names=TRUE,row.names=FALSE)
# selected variables
finaltest_select <- dplyr::select(finaltest, c(game_id, season, team1_id, team2_id,
                                      team1_seed, team2_seed, team1_arate, team2_arate, 
                                      team1_oe, team2_oe, team1_de, team2_de,
                                      team1_pt_school_s16, team2_pt_school_s16, 
                                      team1_pt_overall_s16, team2_pt_overall_s16,
                                      team1_pt_team_season_winrate, team2_pt_team_season_winrate, 
                                      team1_pt_coach_season_winrate, team2_pt_coach_season_winrate))
# write files - with 20 selected variables (no team1WoL, team1WoLI)
write.table(finaltest_select,'/Users/yanqishi/Documents/Data Scientist/March Madness/final/test_data_select.csv',sep=',',col.names=TRUE,row.names=FALSE)

# generate perdiction probability
final_pred <- predict(final_model, newdata=finaltest_select, type='prob')
probability <- final_pred[1:dim(finaltest_select)[1]]
submit <- data.frame(finaltest_select$game_id,probability)
colnames(submit) <- c('game_id','prob')
write.table(submit, '/Users/yanqishi/Documents/Data Scientist/March Madness/final/GUCCIGANG.csv',sep=',',col.names=TRUE,row.names=FALSE)

