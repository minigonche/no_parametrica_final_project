library('depth')
library('MASS')
library('rgl')
library('ggplot2')
library('ddalpha')
library('aplpack')


set.seed(159); 
mu1 <- c(0,0); 
mu2 <- c(6,0); 
sigma <- matrix(c(1,0,0,1), nc = 2)
sigma2 <- matrix(c(8,0,0,3), nc = 2)

mixbivnorm <- rbind(mvrnorm(300, mu1, sigma), mvrnorm(80, mu2, sigma))
singlenorm = mvrnorm(300, mu1, sigma)
uniform = matrix(runif(500, min = -1, max = 1), nc = 2)


#-- COMPARACION PROFUNDIDADES
plot_depth_countour(singlenorm, method = 'Oja')

h = 8
w = 8
#-- VER CONTORNOS

#exp1
pdf("~/Dropbox/Universidad/Materias/Estadistica No parametrica/proyecto_final/exp1_comp.pdf",  height=h, width=w)
plot_depths_with_countours(singlenorm)
dev.off()


pdf("~/Dropbox/Universidad/Materias/Estadistica No parametrica/proyecto_final/exp2_comp.pdf",  height=h, width=w)
plot_depths_with_countours(mixbivnorm, high= 'green')
dev.off()


pdf("~/Dropbox/Universidad/Materias/Estadistica No parametrica/proyecto_final/exp3_comp.pdf",  height=h, width=w)
plot_depths_with_countours(uniform, high = 'red')
dev.off()




#-- CONVEX HULL

convex_hull_peeling(mixbivnorm, hulls = 10)


#-- VER CUANTILES

mean = 0
sd = 1
q1 = 0.1
q2 = 0.9
df = 2

distrub = rnorm(300, mean = mean, sd = sd)
act_l = qnorm(q1,mean = mean, sd = sd)
act_h = qnorm(q2,mean = mean, sd = sd)

distrub = rt(300, df = df)
act_l = qt(q1,df = df)
act_h = qt(q2,df = df)

distrub = rcauchy(1000)
act_l = qcauchy(q1)
act_h = qcauchy(q2)

compare_quantiles(data = distrub , act_l = act_l, act_h = act_h)



rescale = function(data)
{
  return((data - min(data))/(max(data) - min(data)))
}

#Method that receives a bivariate set and plots it with respect to its center
plot_depths = function(data, title = 'Comparacion medidas de profundidad', low="gray50", high="magenta", rescale = TRUE )
{ 
  n = nrow(data)
  t_depth = apply(data, 1, function(x) get_depth(c(x),data, method = 'Tukey') )
  median_t = get_center(data, method = 'Tukey')
  l_depth = apply(data, 1, function(x) get_depth(c(x),data, method = 'Liu') )
  median_l = get_center(data, method = 'Liu')
  o_depth = apply(data, 1, function(x) get_depth(c(x),data, method = 'Oja') )
  median_o = get_center(data, method = 'Oja')
  m_depth = apply(data, 1, function(x) get_depth(c(x),data, method = 'Mahalanobis') )
  median_m = get_center(data, method = 'Mahalanobis')
  
  
  if(rescale)
  {
    t_depth = rescale(t_depth)
    l_depth = rescale(l_depth)
    o_depth = rescale(o_depth)
    m_depth  = rescale(m_depth)
  }
  
  
  d_t = data.frame(x = data[,1], y = data[,2], Type = rep('Tukey') , Depth = t_depth, Median_x = median_t[1], Median_y = median_t[2])
  d_l = data.frame(x = data[,1], y = data[,2], Type = rep('Liu') , Depth = l_depth, Median_x = median_l[1], Median_y = median_l[2])
  d_o = data.frame(x = data[,1], y = data[,2], Type = rep('Oja') , Depth = o_depth, Median_x = median_o[1], Median_y = median_o[2])
  d_m = data.frame(x = data[,1], y = data[,2], Type = rep('Mahalanobis') , Depth = m_depth, Median_x = median_m[1], Median_y = median_m[2])
  
    
    
  dfinal = rbind(d_t,d_l,d_o, d_m)
  
  p = ggplot(dfinal)
  p = p + geom_point( aes(x = x, y = y, size = 4*Depth, colour = Depth))
  p = p + scale_color_gradient(low=low, high=high)
  p = p + geom_point( aes(x = Median_x, y = Median_y), colour = 'black')
  p = p + facet_wrap( ~ Type, scales = "free")
  p = p + ggtitle(title)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

#Method that receives a bivariate set and plots it with respect to its center, and plots countours
plot_depths_with_countours = function(data, title = 'Comparacion medidas de profundidad', low="gray50", high="magenta", rescale = TRUE )
{ 
  n = nrow(data)
  t_depth = apply(data, 1, function(x) get_depth(c(x),data, method = 'Tukey') )
  median_t = get_center(data, method = 'Tukey')
  l_depth = apply(data, 1, function(x) get_depth(c(x),data, method = 'Liu') )
  median_l = get_center(data, method = 'Liu')
  o_depth = apply(data, 1, function(x) get_depth(c(x),data, method = 'Oja') )
  median_o = get_center(data, method = 'Oja')
  m_depth = apply(data, 1, function(x) get_depth(c(x),data, method = 'Mahalanobis') )
  median_m = get_center(data, method = 'Mahalanobis')
  
  
  #countours
  methods = c('Tukey','Liu','Oja','Mahalanobis')
  mins = c(min(t_depth),min(l_depth), min(o_depth), min(m_depth))
  maxs = c(max(t_depth),max(l_depth), max(o_depth), max(m_depth))

  n = 4
  centers_x = c(median_t[1], median_l[1], median_o[1], median_m[1])
  centers_y = c(median_t[2], median_l[2], median_o[2], median_m[2])
  countours = data.frame(x = c(),y = c(), Type = c(), Radius = c())
  
  for(j in 1:length(methods))
  {
    #countours
    method = methods[j]
    min_depth = mins[j]
    max_depth = maxs[j]
    center = c(centers_x[j],centers_y[j])
    
    circles = seq(min_depth, max_depth, (max_depth-min_depth)/(n+1))
    
    for( i in 2:(length(circles)-1))
    {
      radius = circles[i]
      coun = get_contour_circle(data, center = center, method = method, radius = radius, breaks = 100)
      temp = data.frame(x = coun[,1],y = coun[,2], Type = method , Radius = i - 1)
      countours = rbind(countours, temp)
      
    }
  }
  
  if(rescale)
  {
    t_depth = rescale(t_depth)
    l_depth = rescale(l_depth)
    o_depth = rescale(o_depth)
    m_depth  = rescale(m_depth)
  }
  
  d_t = data.frame(x = data[,1], y = data[,2], Type = rep('Tukey') , Depth = t_depth, Median_x = median_t[1], Median_y = median_t[2])
  d_l = data.frame(x = data[,1], y = data[,2], Type = rep('Liu') , Depth = l_depth, Median_x = median_l[1], Median_y = median_l[2])
  d_o = data.frame(x = data[,1], y = data[,2], Type = rep('Oja') , Depth = o_depth, Median_x = median_o[1], Median_y = median_o[2])
  d_m = data.frame(x = data[,1], y = data[,2], Type = rep('Mahalanobis') , Depth = m_depth, Median_x = median_m[1], Median_y = median_m[2])
  
  
  dfinal = rbind(d_t,d_l,d_o, d_m)
  
  p = ggplot(dfinal)
  p = p + geom_point( aes(x = x, y = y, size = 4*Depth, colour = Depth))
  p = p + scale_color_gradient(low=low, high=high)
  p = p + geom_point( aes(x = Median_x, y = Median_y), colour = 'black')
  p = p + geom_path(data = countours[which(countours$Radius ==1),], aes(x = x, y = y ), color = 'blue')
  p = p + geom_path(data = countours[which(countours$Radius ==2),], aes(x = x, y = y ), color = 'blue')
  p = p + geom_path(data = countours[which(countours$Radius ==3),], aes(x = x, y = y ), color = 'blue')
  p = p + geom_path(data = countours[which(countours$Radius ==4),], aes(x = x, y = y ), color = 'blue')
  p = p + facet_wrap( ~ Type, scales = "free")
  p = p + ggtitle(title)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  
  
  return(p)
}

get_center = function(data, method = 'Tukey')
{
  if(method == 'Mahalanobis')
    return(c( mean(data[,1]), mean(data[,2])))
  
  return(med(data, method = method)$median)
  
}

get_depth = function(x,data, method = 'Tukey')
{
  if(method == 'Mahalanobis')
    return(depth.Mahalanobis(x,data))
  
  return(depth(c(x), data, method = method))
}


get_contour_circle = function(data, center = NULL, method = 'Tukey', radius = 0.2, breaks = 10,  eps = 0.0005)
{
 
  if(is.null(center))
  {
    center = get_center(data, method = method)
  }
  center_depth = get_depth(c(center), data, method = method)
  if(radius >= center_depth)
  {
    stop('Impossible length')
  }
  thetas = seq(0, 2*pi, 2*pi/breaks)
  
  #Max length
  max_len = max(apply(data, 1, function(x) sqrt(sum((x - center) ^ 2)) ))
  countour_points = data.frame(x = c(), y = c())
  
  for(theta in thetas)
  {
    #Start binary search
    
    max_len_temp = max_len
    min_len_temp = 0
    len = (max_len_temp-min_len_temp)/2
    point = NULL
    prev_val = Inf
    while(is.null(point))
    {
      temp_point = center + len*c(cos(theta),sin(theta))
      val = get_depth(c(temp_point), data, method = method)
      if(abs(val - radius) < eps ){
        point = temp_point
      }else if(val > radius){
        min_len_temp = len
      }else{
        max_len_temp = len
      }
      
      if(abs(prev_val-val) < eps)
      {
        point = temp_point
      }
      prev_val = val
      len = min_len_temp + (max_len_temp - min_len_temp)/2
      
    }
    
    countour_points = rbind(countour_points,data.frame(x = point[1],y = point[2]))
    
  }
  
  return(countour_points)
  
  
}


plot_depth_countour = function(data, method = 'Tukey', title = '', low="gray50", high="green")
{
  
  depth_data = apply(data, 1, function(x) get_depth(c(x),data, method = method) )
  center = get_center(data, method = method )
  
  min_depth = min(depth_data)
  max_depth = max(depth_data)
  
  n = 4
  circles = seq(min_depth, max_depth, (max_depth-min_depth)/(n+1))
  
  countours = list()
  for( i in 2:(length(circles)-1))
  {
    radius = circles[i]
    coun = get_contour_circle(data, center = center, method = method, radius = radius, breaks = 100)
    countours[i-1] = list(coun)
    
  }
  

  dfinal = data.frame(x = data[,1], y = data[,2] , Depth = depth_data, Median_x = center[1], Median_y = center[2])
  
  p = ggplot(dfinal)
  p = p + geom_point( aes(x = x, y = y, size = 4*Depth, colour = Depth))
  p = p + scale_color_gradient(low=low, high=high)
  #p = p + scale_colour_gradientn(colours = terrain.colors(10))
  p = p + geom_point( aes(x = Median_x, y = Median_y), colour = 'black', size = 3)

  for(coun in countours)
  {
    p = p + geom_path(data = coun, aes(x = x, y = y), colour = 'blue')
  }
  
  p = p + ggtitle(title)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
  
}


compare_quantiles = function(data, act_l, act_h, low_q = 0.1, high_q = 0.9)
{
  j = 1
  samp = data.frame( x = c( act_l, act_h ),
                       y = c(j,j), Type = rep('Actual'))
  
  methods = c('Tukey', 'Oja', 'Liu')
  for(i in 1:length(methods) )
  {
    method = methods[i]
    depths_temp = sapply(data,function(x) get_depth(x,data, method = method))
    q = quantile(depths_temp, probs = c(low_q, high_q))
    sub = data[which(depths_temp >= q[1] & depths_temp <= q[2])]
    x = c(min(sub),max(sub))
    y = c(j+i,j+i)
    
    samp_temp = data.frame( x = x, y = y, Type = rep(method))
    
    samp = rbind(samp, samp_temp)
    
  }
  
  p = ggplot(data = samp)
  p = p + geom_path(aes(x = x, y = y, colour = Type), size = 2)
  p = p + geom_point(aes(x = x, y = rep(0,8), colour = Type, size = 0.7/y))
  p = p + ggtitle('Comparacion cuantiles con medidas de profundidad') + xlab('Posicion Cuantil')
  p = p + theme(plot.title = element_text(hjust = 0.5))
  p
  
}



convex_hull_peeling = function(data, hulls = 10)
{
  plothulls(data, n.hull = hulls)
}



