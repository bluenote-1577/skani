#Models used in the preprint of skani. Trained on only Nayfach et al (2021) GEM data.
MODEL=./skani-mummer-train/v0.1.0-paper-models/0.19525746-3-195-0.06.model
#MODEL=./skani-mummer-train/v0.1.0-paper-models/nayfach_half.model
MODEL_C200=./skani-mummer-train/v0.1.0-paper-models/0.19861849-3-195-0.089999996.model
SRC=src/model.rs

echo $MODEL
echo $MODEL_C200

echo 'pub const MODEL:&str = r#"' > $SRC
cat $MODEL >> $SRC
echo '"#;' >> $SRC
printf '\n' >> $SRC

echo 'pub const MODEL_C200:&str = r#"' >> $SRC
cat $MODEL_C200 >> $SRC
echo '"#;' >> $SRC
printf '\n' >> $SRC
