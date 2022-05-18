function graph_obj = my_rectangle(down_left,up_right)

width  = up_right(1)-down_left(1);
height = up_right(2)-down_left(2);
graph_obj = rectangle('Position',[down_left',width,height]);

end