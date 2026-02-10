# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/b.lind/anaconda3_2025/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/b.lind/anaconda3_2025/etc/profile.d/conda.sh" ]; then
        . "/home/b.lind/anaconda3_2025/etc/profile.d/conda.sh"
    else
        export PATH="/home/b.lind/anaconda3_2025/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

if [ -d "$HOME/pythonimports" ]; then
    export PYTHONPATH="$HOME/pythonimports:$PYTHONPATH"
fi

