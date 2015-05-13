# -*- coding: utf-8 -*-
require 'rake/clean'  # for auto-cleaning
require 'open3'
require 'benchmark'
require_relative 'run_data'
require_relative 'tools/run_iggraph'

namespace :iggraph do 
  num_proc = 8
  
  smab_path = "#{ENV['HOME']}/git_repos/smab_lib/bin/smAb"
  data_f = "data/smab_data_mut10.fa"
  
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

  desc 'run unsupervised on smAb data'
  task :run_smab_data do 
    unless File.exists?("smab_runs_unsup.csv")
      tools = ["igblast", "iggraph"]
      File.open("smab_runs.csv", "w") do |f|
        run_data_unsupervised("./data/smab_data_mut10.fa", f, 
                            "./smab_run_pred", tools)
      end
    end
  end

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
  
  desc 'run on stanford_s22 data'
  task :test_stanford_s22_data do 
    unless File.exists?("stanford_s22_runs.csv")
      #tools = ["igblast", "iggraph", "ihmmune"]
      tools = ["igblast", "iggraph"]
      #tools = ["igblast"]
      File.open("stanford_s22_runs.csv", "w") do |f|
        run_data_unsupervised("./data/Stanford_S22_upcase.fasta",
                              f, 
                              "./stanford_s22_run_pred",
                              tools, 
                              "allele")
      end
    end
  end
  #CLEAN << "stanford_s22_runs.csv"
  #CLEAN << Dir.glob("stanford_s22_run_pred*")
  
  
  desc 'compare Stanford_S22 data predictions'
  task :cmp_stanford_s22_data => :test_stanford_s22_data do
    out_file_base = "stanford_s22_data_comparison"
    type = "allele"
    unless !Dir.glob(out_file_base+"*").empty?
      #
      puts 'hi'
      pred_lst = Dir.glob("./stanford_s22_run_pred*")
      compare_unsupervised(pred_lst, out_file_base, "allele")
    end
  end
  CLOBBER << Dir.glob("stanford_s22_data_comparison*")

  desc 'compare all Stanford_S22 tools'  
  s22_out_file_base = "stanford_s22_all_comparison"
  task :cmp_all_stanford_s22 do     
    name_re = Regexp.new(/(\w+)\_S22\_results\.txt$/)
    unless !Dir.glob(s22_out_file_base + "*").empty?
      ["gene", "allele"].each do |y|
        out_file_base_s = [s22_out_file_base, y].join("_")
        pred_lst = Dir.glob("./data/s22_results/*results.txt")        
        #puts out_file_base_s
        #
        name_lst = pred_lst.map{|s| name_re.match(s)[1]}
        pred_lst.push("./data/stanford_s22_run_pred_iggraph.tab")
        name_lst.push("IgGraph")
        # puts pred_lst.to_s
        # puts name_lst.to_s
        compare_unsupervised(pred_lst, 
                             out_file_base_s, 
                             y, 
                             name_lst)
      end
    end
  end
  CLEAN << Dir.glob(s22_out_file_base + "*")

  desc 'plot all Stanford_S22'
  task :plot_all_stanford_s22 => :cmp_all_stanford_s22 do 
    unless !Dir.glob("stanford_s22_mat_*.pdf").empty?
      ["V", "D", "J", "total"].each do |seg|
        cmd = "./plot_matrix.R stanford_s22_all_comparison_allele_#{seg}.csv stanford_s22_mat_#{seg}.pdf half_mat"
        Open3.capture3(cmd)
      end
    end
  end
  CLEAN << Dir.glob("stanford_s22_mat_*.pdf")

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

  desc 'label mouse Ig-seq'
  task :run_mouse_igseq do
    unless !Dir.glob("mouse_label*.tab").empty?
      mouse_igseq = "../data/mouse_igseq.fa"
      cmd = "./label_ig.rb #{mouse_igseq} -p 7 -o ../mouse_label -g mouse"
      out,err,pip = Open3.capture3(cmd, :chdir => "tools")
      puts out
      File.open("mouse_label_runtime.txt","w"){|f| f.puts out}
    end
  end
  
  desc 'label human Ig-seq'
  task :run_human_igseq do
    unless !Dir.glob("human_label*.tab").empty?
      human_igseq = 
        "#{ENV['HOME']}/data/ig_seq/gen_data/7_SAM15574987_HС_naive.clusters.fa"
      cmd = "./label_ig.rb #{human_igseq} -p 7 -o ../human_label -g human"
      out,err,pip = Open3.capture3(cmd, :chdir => "tools")
      puts out
      File.open("human_label_runtime.txt","w"){|f| f.puts out}
    end
  end
  
end
