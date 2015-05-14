# -*- coding: utf-8 -*-
require 'rake/clean'  # for auto-cleaning
require 'open3'
require 'benchmark'
require_relative 'run_data'
require_relative 'tools/run_iggraph'

namespace :iggraph do 
  #
  # unsupervised tasks
  #
  load 'unsupervised.rake'
  #
  #
  num_proc = 8
  #
  #
  smab_path = "#{ENV['HOME']}/git_repos/smab_lib/bin/smAb"
  data_f = "data/smab_data_mut10.fa"
  #
  #
  desc 'grid search over VJ/D k-mers for alleles'
  task :grid_search_k_allele do
    unless File.exists?("param_out_alleles.csv")
      cmd = "./param_search.rb #{data_f} param_out_alleles.csv allele #{num_proc}"
      Open3.capture3(cmd)      
    end
  end
  #CLOBBER << "param_out_alleles.csv" # clobber calls clean

  desc 'plot VJ/D k-mer performance for alleles'
  task :plot_k_allele => :grid_search_k_allele do
    unless File.exists?("mat_allele.pdf")      
      check_cmd = "#{smab_path} check -f #{data_f}"
      out_r,err,pip = Open3.capture3(check_cmd)
      max_rate = out_r.chomp.to_f
      cmd = "./plot_matrix.R param_out_alleles.csv mat_allele.pdf heatmap #{max_rate}"
      Open3.capture3(cmd)
    end
  end
  CLEAN << "mat_allele.pdf"
  
  desc 'grid search over VJ/D k-mers for gene'
  task :grid_search_k_gene do
    unless File.exists?("param_out_gene.csv")
      cmd = "./param_search.rb #{data_f} param_out_gene.csv gene #{num_proc}"
      Open3.capture3(cmd)
    end
  end
  #CLOBBER << "param_out_gene.csv"

  desc 'plot VJ/D k-mer performance for genes'
  task :plot_k_gene => :grid_search_k_gene do
    unless File.exists?("mat_gene.pdf")
      check_cmd = "#{smab_path} check -f #{data_f}"
      out_r,err,pip = Open3.capture3(check_cmd)
      max_rate = out_r.chomp.to_f
      cmd = "./plot_matrix.R param_out_gene.csv mat_gene.pdf heatmap #{max_rate}"
      Open3.capture3(cmd)
    end
  end
  CLEAN << "mat_gene.pdf"

  desc 'run on smAb data'
  task :test_smab_data do 
    unless File.exists?("smab_runs.csv")
      tools = ["igblast", "iggraph", "ihmmune"]
      File.open("smab_runs.csv", "w") do |f|
        run_data_supervised("./data/smab_data_mut10.fa", f, 
                            "./smab_run_pred", tools, "allele")
        #run_data_supervised("../tests/data/ten_ab.fa", f, "./smab_run_pred", "allele")
      end
    end
  end
  #CLEAN << "smab_runs.csv"
  #CLEAN << Dir.glob("smab_run_pred*")
  
  desc 'compare smAb data predictions'
  task :cmp_smab_data => :test_smab_data do 
    out_file_base = "smab_data_comparison"
    type = "allele"
    unless !Dir.glob(out_file_base+"*").empty?
      #
      pred_lst = Dir.glob("./smab_run_pred*")
      compare_unsupervised(pred_lst, out_file_base, "allele")      
    end
  end
  #CLOBBER << Dir.glob("smab_data_comparison*")
  
  desc 'benchmark iggraph options'
  task :benchmark_iggraph_options do 
    Benchmark.bm() do |bm|
      fasta_f = "data/smab_data_mut10.fa"
      fasta_p = File.expand_path(File.dirname(__FILE__)) + "/" + fasta_f      
      bm.report("std") {
        rigg = RunIgGraph.new(fasta_f)
        rigg.compute
        rigg.write_preds("/tmp/test_runigg_std.tab")
        rigg.cleanup
      }
      bm.report("std+no_cdr3") {
        rigg = RunIgGraph.new(fasta_f)
        rigg.no_cdr3 = true
        rigg.compute
        rigg.write_preds("/tmp/test_runigg_nocdr3.tab")
        rigg.cleanup
      }
      # bm.report("prob") {
      # }
    end
  end
  
end
