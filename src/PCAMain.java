import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.mllib.linalg.Matrix;
import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;
import org.apache.spark.mllib.linalg.distributed.RowMatrix;

import java.util.Arrays;
import java.util.List;

/**
 * Created by 李凌耀 on 2017/8/4.
 */

/*
 *  PCA算法实现
 *  设初始数据为m条n维矩阵
 *  1、组成n行m列矩阵X
 *  2、零均值化
 *  3、求出协方差矩阵 C=1/m*X*X(T)   （T）表示转置
 *  4、求出协方差矩阵的特征值和特征向量
 *  5、按特征值从大到小从上到下排列矩阵，取前k行组成P
 *  6、Y=PX 即为降维到k维后的数据
 *  输入格式
 *  a 1 2 3 4
 *  b 2 3 4 5
 */
public class PCAMain {
    public static void main(String[] args) {
        SparkConf conf = new SparkConf().setAppName("PCA");
        JavaSparkContext sc = new JavaSparkContext(conf);
        JavaRDD<String> input = sc.textFile("file");

        JavaRDD<Vector> rows = input.map(
                new Function<String, Vector>() {
                    @Override
                    public Vector call(String s) throws Exception {
                        String[] split = s.split(" ");
                        double[] a = new double[split.length-1];
                        for (int i=0;i<a.length;i++)
                        {
                            a[i] = Double.parseDouble(split[i+1]);
                        }
                        return Vectors.dense(a);
                    }
                }
        );

        RowMatrix mat = new RowMatrix(rows.rdd());

        Matrix pc = mat.computePrincipalComponents(2);

        RowMatrix projected = mat.multiply(pc);

        Vector[] collectPartitions = (Vector[])projected.rows().collect();

        List<Vector> data = Arrays.asList(collectPartitions);;

        JavaRDD<Vector> vec = sc.parallelize(data);

        vec.saveAsTextFile("savePath");


    }
}
