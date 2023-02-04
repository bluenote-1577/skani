#MODEL=skani_mummer_train/good_models/0.21994846-4-180-0.08.model
MODEL=skani-mummer-train/good_models/0.19134241-3-175-0.12.model
MODEL_C200=skani-mummer-train/good_models/0.18635842-7-230-0.08.model
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
