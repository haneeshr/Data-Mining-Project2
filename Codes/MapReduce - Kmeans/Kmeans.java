import java.io.*;
import java.util.StringTokenizer;
import java.util.*;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.*;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;

public class Kmeans {

  static int k=5;
  static double[][] centroids;
  static double[][] temp;
  static String filnename = "centroids.txt";
  static int dimension;

  public static class TokenizerMapper
       extends Mapper<Object, Text, IntWritable, Text>{

    private final static IntWritable one = new IntWritable(1);

    public void map(Object key, Text value, Context context
                    ) throws IOException, InterruptedException {
      
      for(String text : value.toString().split("\n")){
        String[] split = text.toString().split("\t");
        double[] point = new double[split.length-2];

        for(int i=0; i<point.length;i++){
          point[i] = Double.parseDouble(split[i+2]);
        }

        int closest = 0;
        double min = Double.MAX_VALUE;

        for(int i=0; i<centroids.length; i++){
          double dist = 0.0;
          for(int j=0; j<point.length; j++){
            dist = dist + (centroids[i][j] - point[j])*(centroids[i][j] - point[j]);
          }
          dist = Math.pow(dist, 0.5);
          if(dist < min){
            min = dist;
            closest = i;
          }
        }
        IntWritable result = new IntWritable(closest);
        context.write(result, value);
      }
      
    }
  }

  public static class IntSumReducer
       extends Reducer<IntWritable,Text,IntWritable,Text> {

    public void reduce(IntWritable key, Iterable<Text> values,
                       Context context
                       ) throws IOException, InterruptedException {

      double[] sum = new double[dimension];


      double count=0;
      for(Text text : values){
        String[] split = text.toString().split("\t");
        for(int i=0;i<split.length-2;i++){
          sum[i] += Double.parseDouble(split[i+2]);
        }
        count++;
        context.write(new IntWritable(Integer.parseInt(split[0])), new Text(key.get()+""));
      }

      String res = "";
      for(int i=0; i<sum.length; i++){
        sum[i] = sum[i]/count;
        if(i!=0)res+="\t";
        res+=sum[i];
        temp[key.get()][i] = sum[i];
      }
    }
  }

  public static void main(String[] args) throws Exception {

    int count=-1;


    k=Integer.parseInt(args[2]);
    count=Integer.parseInt(args[3]);

    Configuration conf;
    Job job;
    BufferedReader bufferedReader = new BufferedReader(new FileReader(filnename));
    String sCurrentLine = bufferedReader.readLine();
    dimension = sCurrentLine.split("\t").length;


    centroids = new double[k][dimension];
    temp = new double[k][dimension];

    for(int i=0; i<k; i++){
      String split[] = sCurrentLine.split("\t");
      for(int j=0; j<dimension; j++){
        centroids[i][j] = Double.parseDouble(split[j]);
      }
      sCurrentLine = bufferedReader.readLine();
    }

    

    do{
      conf = new Configuration();
      job = Job.getInstance(conf, "word count");
      job.setJarByClass(Kmeans.class);
      job.setMapperClass(TokenizerMapper.class);
      // job.setCombinerClass(IntSumReducer.class);
      job.setReducerClass(IntSumReducer.class);
      job.setOutputKeyClass(IntWritable.class);
      job.setOutputValueClass(Text.class);
      FileInputFormat.addInputPath(job, new Path(args[0]));
      FileOutputFormat.setOutputPath(job, new Path(args[1]));
      job.waitForCompletion(true);

      if(converged(temp))break;
      if(count!=1)FileSystem.get(conf).delete(new Path(args[1]), true);
      for(int i = 0; i < temp.length; i++){
        for(int j = 0; j < temp[0].length ; j++){
            centroids[i][j] = temp[i][j];
        }
      }
      if(count>0)count--;
    }while(true && count!=0);
  }

  public static boolean converged(double[][] temp){

    for(int i=0; i<temp.length; i++){
      for(int j=0; j<temp[0].length; j++){
        if(temp[i][j]!=centroids[i][j]){
          return false;
        }
      }
    }
    return true;
  }
}