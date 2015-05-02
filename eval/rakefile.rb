require 'rake/clean'  # for auto-cleaning
require 'open3'

namespace :iggraph do 
  num_proc = 8
  
  smab_path = "#{ENV['HOME']}/git_repos/smab_lib/bin/smAb"
  data_f = "data/smab_data_mut10.fa"
  
  desc 'grid search over VJ/D k-mers for alleles'
  task :grid_search_k_allele do
    unless File.exists?("param_out_alleles.csv")
      #data_f = "data/smab_data_mut10.fa"
      cmd = "./param_search.rb #{data_f} param_out_alleles.csv allele #{num_proc}"
      Open3.capture3(cmd)      
    end
  end
  CLOBBER << "param_out_alleles.csv" # clobber calls clean

  desc 'plot VJ/D k-mer performance for alleles'
  task :plot_k_allele => :grid_search_k_allele do
    unless File.exists?("mat_allele.pdf")      
      check_cmd = "#{smab_path} check -f #{data_f}"
      out_r,err,pip = Open3.capture3(check_cmd)
      max_rate = out_r.chomp.to_f
      #puts "RATE = #{max_rate}"
      cmd = "./plot_matrix.R param_out_alleles.csv mat_allele.pdf heatmap #{max_rate}"
      #puts cmd
      Open3.capture3(cmd)
    end
  end
  CLEAN << "mat_allele.pdf"
  
  desc 'grid search over VJ/D k-mers for gene'
  task :grid_search_k_gene do
    unless File.exists?("param_out_gene.csv")
      #data_f = "data/smab_data_mut10.fa"
      cmd = "./param_search.rb #{data_f} param_out_gene.csv gene #{num_proc}"
      Open3.capture3(cmd)
    end
  end
  CLOBBER << "param_out_gene.csv"

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
  
end
