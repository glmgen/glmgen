
cd ~/tmp_tf/
/bin/cp -f tf_admm.c tf_admm_glm.c tf.h tf_maxlam.c ~/code/glmgen/glmgen/R_pkg/glmgen/src/tf/
/bin/cp -f utils.h diag_times_sparse.c ~/code/glmgen/glmgen/R_pkg/glmgen/src/utils/
cd -

# git checkout glmgen/src/tf/tf.h glmgen/src/tf/tf_admm.c glmgen/src/tf/tf_admm_glm.c glmgen/src/tf/tf_maxlam.c glmgen/src/utils/utils.h
