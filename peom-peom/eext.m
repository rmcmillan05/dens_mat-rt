function y = eext(t)

global E0
global start_delta_id
global step_id

if step_id == start_delta_id
    y = E0;
else
    y = 0;
end